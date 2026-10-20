#' Import a TCGA/GDC project
#'
#' @description
#' Downloads and assembles a Genomic Data Commons project into a single
#' \linkS4class{OncoSeqExperiment}: STAR gene-expression counts, the
#' indexed clinical table with harmonised overall-survival columns, and
#' optionally simple somatic mutations and gene-level copy-number scores.
#' Layers are keyed on the TCGA patient barcode so that expression,
#' genomic and clinical data for one tumour line up.
#'
#' @details
#' The GDC retired the `"HTSeq - Counts"` workflow when it moved to Data
#' Release 32; this function requests `"STAR - Counts"`, which is the
#' current gene-expression workflow. Downloaded files are cached by
#' `TCGAbiolinks` in `directory`, so repeating a call is inexpensive.
#'
#' Sample barcodes are truncated to the 12-character patient barcode so
#' that expression, mutation, copy-number and clinical layers share
#' identifiers. When a patient contributes more than one aliquot the
#' first is retained and a message reports how many were dropped.
#'
#' @param project GDC project identifier, for example `"TCGA-BLCA"` for
#'   bladder urothelial carcinoma.
#' @param sampleType Character vector of sample types to keep, for example
#'   `c("Primary Tumor", "Solid Tissue Normal")`. `NULL` keeps all.
#' @param mutations Logical; also download simple somatic mutations.
#' @param cnv Logical; also download gene-level copy-number scores.
#' @param clinical Logical; attach the indexed clinical table.
#' @param directory Download cache directory passed to
#'   [TCGAbiolinks::GDCdownload()].
#' @param files.per.chunk Number of files per download chunk; lower this on
#'   unstable connections.
#'
#' @return An \linkS4class{OncoSeqExperiment} with a `counts` assay (plus
#'   `tpm_unstrand` and `fpkm_unstrand` when the GDC supplies them),
#'   clinical columns in `colData`, and the mutation and copy-number layers
#'   when requested.
#'
#' @section Requirements:
#' Needs the `TCGAbiolinks` package and network access to
#' `api.gdc.cancer.gov`.
#'
#' @examples
#' # Downloads several hundred megabytes; not run during package checks.
#' \dontrun{
#' blca <- importGDC("TCGA-BLCA", sampleType = "Primary Tumor",
#'                   mutations = TRUE, cnv = TRUE)
#' blca
#' }
#'
#' @seealso [importCBioPortal()] for a lighter, API-based alternative.
#' @export
importGDC <- function(project = "TCGA-BLCA", sampleType = NULL,
                      mutations = FALSE, cnv = FALSE, clinical = TRUE,
                      directory = "GDCdata", files.per.chunk = 50L) {
    .requireNamespaces("TCGAbiolinks", "importGDC()")

    query <- TCGAbiolinks::GDCquery(
        project = project,
        data.category = "Transcriptome Profiling",
        data.type = "Gene Expression Quantification",
        workflow.type = "STAR - Counts",
        sample.type = sampleType)
    TCGAbiolinks::GDCdownload(query, directory = directory,
                              files.per.chunk = files.per.chunk)
    se <- TCGAbiolinks::GDCprepare(query, directory = directory)

    assayList <- list(counts = SummarizedExperiment::assay(se, "unstranded"))
    for (extra in c("tpm_unstrand", "fpkm_unstrand")) {
        if (extra %in% SummarizedExperiment::assayNames(se)) {
            assayList[[extra]] <- SummarizedExperiment::assay(se, extra)
        }
    }

    cd <- as.data.frame(SummarizedExperiment::colData(se))
    cd$sample_barcode <- colnames(se)
    cd$patient <- .tcgaPatientBarcode(colnames(se))

    dedup <- .dedupByPatient(cd$patient)
    if (sum(!dedup) > 0L) {
        message("dropping ", sum(!dedup),
                " additional aliquot(s) so that one sample remains per patient")
    }
    assayList <- lapply(assayList, function(m) m[, dedup, drop = FALSE])
    cd <- cd[dedup, , drop = FALSE]
    rownames(cd) <- cd$patient
    assayList <- lapply(assayList, function(m) {
        colnames(m) <- cd$patient
        m
    })

    rd <- as.data.frame(SummarizedExperiment::rowData(se))

    if (clinical) {
        clin <- try(TCGAbiolinks::GDCquery_clinic(project, type = "clinical"),
                    silent = TRUE)
        if (!inherits(clin, "try-error")) {
            clin <- .harmoniseSurvivalColumns(clin, idColumn = "submitter_id")
            idx <- match(rownames(cd), clin$submitter_id)
            keepCols <- setdiff(colnames(clin), colnames(cd))
            cd <- cbind(cd, clin[idx, keepCols, drop = FALSE])
            rownames(cd) <- colnames(assayList[[1L]])
        } else {
            warning("clinical data could not be retrieved for ", project,
                    call. = FALSE)
        }
    }

    mutDF <- NULL
    if (mutations) {
        mutDF <- .gdcMutations(project, directory = directory,
                               keep = rownames(cd))
    }

    cnvMat <- NULL
    if (cnv) {
        cnvMat <- .gdcCNV(project, directory = directory, keep = rownames(cd))
    }

    OncoSeqExperiment(
        assays = assayList,
        colData = cd,
        rowData = rd,
        mutations = mutDF,
        cnv = cnvMat,
        provenance = list(source = "GDC", project = project,
                          workflow = "STAR - Counts",
                          retrieved = Sys.time()))
}

.gdcMutations <- function(project, directory, keep) {
    q <- try(TCGAbiolinks::GDCquery(
        project = project,
        data.category = "Simple Nucleotide Variation",
        data.type = "Masked Somatic Mutation",
        access = "open"), silent = TRUE)
    if (inherits(q, "try-error")) {
        warning("no open-access mutation data for ", project, call. = FALSE)
        return(NULL)
    }
    TCGAbiolinks::GDCdownload(q, directory = directory)
    maf <- TCGAbiolinks::GDCprepare(q, directory = directory)
    maf <- as.data.frame(maf)
    maf$sample <- .tcgaPatientBarcode(maf$Tumor_Sample_Barcode)
    maf <- maf[maf$sample %in% keep, , drop = FALSE]
    .normaliseMutationTable(maf)
}

.gdcCNV <- function(project, directory, keep) {
    q <- try(TCGAbiolinks::GDCquery(
        project = project,
        data.category = "Copy Number Variation",
        data.type = "Gene Level Copy Number Scores"), silent = TRUE)
    if (inherits(q, "try-error")) {
        warning("no gene-level copy-number scores for ", project,
                call. = FALSE)
        return(NULL)
    }
    TCGAbiolinks::GDCdownload(q, directory = directory)
    scores <- TCGAbiolinks::GDCprepare(q, directory = directory)
    scores <- as.data.frame(scores)
    geneCol <- intersect(c("Gene Symbol", "gene_name", "Gene.Symbol"),
                         colnames(scores))[1L]
    if (is.na(geneCol)) {
        warning("unexpected copy-number table layout; skipping CNV layer",
                call. = FALSE)
        return(NULL)
    }
    sampleCols <- setdiff(colnames(scores),
        c(geneCol, "Gene ID", "Cytoband", "gene_id", "Gene.ID"))
    mat <- as.matrix(scores[, sampleCols, drop = FALSE])
    rownames(mat) <- make.unique(as.character(scores[[geneCol]]))
    colnames(mat) <- .tcgaPatientBarcode(colnames(mat))
    mat <- mat[, colnames(mat) %in% keep, drop = FALSE]
    mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]
    if (ncol(mat) == 0L) return(NULL)
    mat
}

## TCGA barcodes: TCGA-XX-YYYY-... ; the first three fields identify the patient.
.tcgaPatientBarcode <- function(x) {
    x <- as.character(x)
    vapply(strsplit(x, "-", fixed = TRUE), function(p) {
        if (length(p) >= 3L) paste(p[seq_len(3L)], collapse = "-") else p[1L]
    }, character(1))
}

.dedupByPatient <- function(patient) {
    !duplicated(patient)
}

## Build 'os_time' (days) and 'os_event' (1 = death) from GDC clinical fields.
.harmoniseSurvivalColumns <- function(clin, idColumn = "submitter_id") {
    clin <- as.data.frame(clin)
    status <- if ("vital_status" %in% colnames(clin)) {
        tolower(as.character(clin$vital_status))
    } else rep(NA_character_, nrow(clin))
    event <- ifelse(status == "dead", 1L, ifelse(status == "alive", 0L, NA_integer_))

    death <- .numericColumn(clin, c("days_to_death"))
    follow <- .numericColumn(clin,
        c("days_to_last_follow_up", "days_to_last_followup"))
    time <- ifelse(!is.na(death), death, follow)

    clin$os_time <- time
    clin$os_event <- event
    clin
}

.numericColumn <- function(df, candidates) {
    hit <- intersect(candidates, colnames(df))
    if (length(hit) == 0L) return(rep(NA_real_, nrow(df)))
    suppressWarnings(as.numeric(as.character(df[[hit[1L]]])))
}
