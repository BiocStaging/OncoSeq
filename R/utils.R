#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom methods is as new
#' @importFrom utils head read.delim packageVersion
NULL

## Coerce a data.frame / DataFrame / matrix to DataFrame, keeping row names.
.asDataFrame <- function(x) {
    if (is.null(x)) return(NULL)
    if (is(x, "DataFrame")) return(x)
    DataFrame(x, check.names = FALSE)
}

.showLine <- function(label, values) {
    if (length(values) == 0L) {
        return(paste0(label, "(0):\n"))
    }
    shown <- head(values, 6L)
    extra <- if (length(values) > 6L) " ..." else ""
    paste0(label, "(", length(values), "): ",
           paste(shown, collapse = ", "), extra, "\n")
}

## Column synonyms accepted when a MAF-style table is supplied.
.MUTATION_SYNONYMS <- list(
    sample = c("sample", "Tumor_Sample_Barcode", "sample_id", "sampleId",
               "sampleID", "Sample_ID", "patientId"),
    gene = c("gene", "Hugo_Symbol", "hugoGeneSymbol", "gene_symbol",
             "Gene", "symbol"),
    variant_class = c("variant_class", "Variant_Classification",
                      "mutationType", "variantClassification",
                      "consequence"),
    variant_type = c("variant_type", "Variant_Type", "variantType"),
    chromosome = c("chromosome", "Chromosome", "chr"),
    start = c("start", "Start_Position", "startPosition"),
    end = c("end", "End_Position", "endPosition"),
    ref = c("ref", "Reference_Allele", "referenceAllele"),
    alt = c("alt", "Tumor_Seq_Allele2", "variantAllele"),
    protein_change = c("protein_change", "HGVSp_Short", "proteinChange")
)

## Rename recognised columns of a variant table to the OncoSeq vocabulary
## without discarding any of the original columns.
.harmoniseMutationColumns <- function(df) {
    nms <- colnames(df)
    for (target in names(.MUTATION_SYNONYMS)) {
        if (target %in% nms) next
        hit <- intersect(.MUTATION_SYNONYMS[[target]], nms)
        if (length(hit) > 0L) {
            df[[target]] <- df[[hit[1L]]]
        }
    }
    df
}

.normaliseMutationTable <- function(mutations) {
    if (is.null(mutations)) return(DataFrame())
    if (is(mutations, "DataFrame")) {
        df <- mutations
    } else if (is.data.frame(mutations)) {
        df <- DataFrame(mutations, check.names = FALSE)
    } else {
        stop("'mutations' must be a data.frame or DataFrame", call. = FALSE)
    }
    if (nrow(df) == 0L) return(DataFrame())
    df <- .harmoniseMutationColumns(df)
    missing_cols <- setdiff(c("sample", "gene"), colnames(df))
    if (length(missing_cols) > 0L) {
        stop("the mutation table needs column(s) ",
             paste(missing_cols, collapse = " and "),
             "; none of the recognised synonyms were found. Supplied ",
             "columns: ", paste(colnames(df), collapse = ", "),
             call. = FALSE)
    }
    df$sample <- as.character(df$sample)
    df$gene <- as.character(df$gene)
    if ("variant_class" %in% colnames(df)) {
        df$variant_class <- as.character(df$variant_class)
    }
    rownames(df) <- NULL
    df
}

.asCNVSummarizedExperiment <- function(cnv) {
    if (is.null(cnv)) return(NULL)
    if (is(cnv, "SummarizedExperiment")) return(cnv)
    if (is.data.frame(cnv)) cnv <- as.matrix(cnv)
    if (!is.matrix(cnv)) {
        stop("'cnv' must be a matrix, data.frame or SummarizedExperiment",
             call. = FALSE)
    }
    if (is.null(rownames(cnv)) || is.null(colnames(cnv))) {
        stop("the copy-number matrix must have both row and column names",
             call. = FALSE)
    }
    SummarizedExperiment::SummarizedExperiment(
        assays = list(cn = cnv))
}

## Stop with an actionable message when an optional dependency is absent.
.requireNamespaces <- function(pkgs, what) {
    missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1),
                                 quietly = TRUE)]
    if (length(missing_pkgs) > 0L) {
        stop(what, " requires the package(s) ",
             paste(sQuote(missing_pkgs), collapse = ", "),
             ", which are not installed. Install them with ",
             "BiocManager::install(c(",
             paste(paste0('"', missing_pkgs, '"'), collapse = ", "), "))",
             call. = FALSE)
    }
    invisible(TRUE)
}

## Append a step to the provenance log.
.addProvenance <- function(object, step, ...) {
    entry <- list(step = step, time = Sys.time(), params = list(...))
    object@provenance$steps <- c(object@provenance$steps, list(entry))
    object
}

## Store an analysis result in the registry.
.addAnalysis <- function(object, name, result, analysis = name, ...) {
    attr(result, "analysis") <- analysis
    attr(result, "params") <- list(...)
    current <- object@analyses
    current[[name]] <- result
    object@analyses <- current
    object
}

#' Retrieve a stored analysis result
#'
#' @description Extract one result from the analysis registry of an
#'   \linkS4class{OncoSeqExperiment}, with an informative error when the name
#'   is not present.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param name Name of the stored result. When missing, the available names
#'   are returned.
#'
#' @return The stored result object, or a character vector of available names
#'   when `name` is missing.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' getAnalysis(ose)
#' ose <- diffExpress(ose, ~ condition,
#'                    contrast = c("condition", "tumour", "normal"),
#'                    engine = "edgeR", verbose = FALSE)
#' head(getAnalysis(ose, "de"))
#'
#' @export
getAnalysis <- function(object, name) {
    stopifnot(is(object, "OncoSeqExperiment"))
    available <- names(object@analyses)
    if (is.null(available)) available <- character(0)
    if (missing(name)) return(available)
    if (!name %in% available) {
        stop("no analysis named ", sQuote(name), " in this object. ",
             if (length(available) == 0L) "The registry is empty."
             else paste0("Available: ", paste(available, collapse = ", ")),
             call. = FALSE)
    }
    object@analyses[[name]]
}

## Resolve the assay to work on, preferring an explicit name.
.getAssay <- function(object, assay = NULL, default = "counts") {
    nms <- SummarizedExperiment::assayNames(object)
    if (is.null(assay)) {
        assay <- if (default %in% nms) default else nms[1L]
    }
    if (!assay %in% nms) {
        stop("assay ", sQuote(assay), " not found; available: ",
             paste(nms, collapse = ", "), call. = FALSE)
    }
    SummarizedExperiment::assay(object, assay)
}

## Benjamini-Hochberg helper that tolerates all-NA input.
.bh <- function(p) {
    out <- rep(NA_real_, length(p))
    ok <- !is.na(p)
    if (any(ok)) out[ok] <- stats::p.adjust(p[ok], method = "BH")
    out
}
