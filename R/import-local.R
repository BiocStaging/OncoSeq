#' Import counts and sample annotation from local data
#'
#' @description
#' Builds an \linkS4class{OncoSeqExperiment} from a gene-by-sample count
#' matrix and a table of sample annotation held in memory or on disk.
#' Optional somatic-variant and copy-number layers are attached and checked
#' against the expression sample identifiers.
#'
#' @param counts A numeric matrix (genes in rows, samples in columns), a
#'   `data.frame`, or a path to a delimited text file whose first column holds
#'   feature identifiers.
#' @param colData A `data.frame`/\link[S4Vectors:DataFrame-class]{DataFrame} of sample annotation,
#'   or a path to a delimited file whose first column holds sample
#'   identifiers. Rows are matched to the columns of `counts` by name when
#'   names are available, otherwise by position.
#' @param rowData Optional feature annotation, as an object or a file path.
#' @param mutations Optional variant table, as an object or a path to a MAF
#'   or delimited file (see [importMAF()]).
#' @param cnv Optional gene-level copy-number matrix, as an object or a file
#'   path.
#' @param sep Field separator used when reading files.
#' @param intersectSamples When `TRUE` (the default) the object is restricted
#'   to samples present in both `counts` and `colData`. When `FALSE` a
#'   mismatch is an error.
#'
#' @return An \linkS4class{OncoSeqExperiment}.
#'
#' @examples
#' counts <- matrix(rpois(300, 40), nrow = 30,
#'                  dimnames = list(paste0("G", 1:30), paste0("S", 1:10)))
#' meta <- data.frame(condition = rep(c("tumour", "normal"), each = 5),
#'                    row.names = colnames(counts))
#' ose <- importCounts(counts, meta)
#' ose
#'
#' @seealso [importGDC()], [importCBioPortal()], [importTximport()]
#' @export
importCounts <- function(counts, colData = NULL, rowData = NULL,
                         mutations = NULL, cnv = NULL, sep = "\t",
                         intersectSamples = TRUE) {
    counts <- .readIfPath(counts, sep = sep, rowNames = TRUE)
    if (is.data.frame(counts)) counts <- as.matrix(counts)
    if (!is.matrix(counts)) {
        stop("'counts' must be a matrix, data.frame or a path to one",
             call. = FALSE)
    }
    if (!is.numeric(counts)) {
        stop("'counts' must be numeric; non-numeric columns were found. ",
             "Check that feature identifiers are in the row names rather ",
             "than in a data column.", call. = FALSE)
    }
    if (is.null(rownames(counts))) {
        stop("'counts' must have row names giving feature identifiers",
             call. = FALSE)
    }
    if (is.null(colnames(counts))) {
        colnames(counts) <- paste0("sample", seq_len(ncol(counts)))
    }

    cd <- .readIfPath(colData, sep = sep, rowNames = TRUE)
    if (!is.null(cd)) {
        cd <- as.data.frame(cd)
        if (!is.null(rownames(cd)) &&
            length(intersect(rownames(cd), colnames(counts))) > 0L) {
            common <- intersect(colnames(counts), rownames(cd))
            if (!intersectSamples &&
                (length(common) != ncol(counts) || length(common) != nrow(cd))) {
                stop("sample identifiers in 'counts' and 'colData' differ ",
                     "and 'intersectSamples' is FALSE", call. = FALSE)
            }
            if (length(common) == 0L) {
                stop("no sample identifiers shared by 'counts' and 'colData'",
                     call. = FALSE)
            }
            counts <- counts[, common, drop = FALSE]
            cd <- cd[common, , drop = FALSE]
        } else if (nrow(cd) == ncol(counts)) {
            rownames(cd) <- colnames(counts)
        } else {
            stop("'colData' has ", nrow(cd), " rows but 'counts' has ",
                 ncol(counts), " columns, and their names do not match",
                 call. = FALSE)
        }
    }

    rd <- .readIfPath(rowData, sep = sep, rowNames = TRUE)
    if (!is.null(rd)) {
        rd <- as.data.frame(rd)
        if (!is.null(rownames(rd))) {
            common <- intersect(rownames(counts), rownames(rd))
            if (length(common) > 0L) {
                rd <- rd[match(rownames(counts), rownames(rd)), , drop = FALSE]
                rownames(rd) <- rownames(counts)
            }
        }
    }

    mut <- if (is.character(mutations) && length(mutations) == 1L) {
        .readMAFTable(mutations, sep = sep)
    } else mutations

    cnvMat <- .readIfPath(cnv, sep = sep, rowNames = TRUE)
    if (is.data.frame(cnvMat)) cnvMat <- as.matrix(cnvMat)
    if (!is.null(cnvMat)) {
        keep <- intersect(colnames(cnvMat), colnames(counts))
        if (length(keep) == 0L) {
            stop("the copy-number data shares no sample identifiers with ",
                 "the expression data", call. = FALSE)
        }
        cnvMat <- cnvMat[, keep, drop = FALSE]
    }

    if (!is.null(mut)) {
        mutDF <- .normaliseMutationTable(mut)
        if (nrow(mutDF) > 0L) {
            keep <- as.character(mutDF$sample) %in% colnames(counts)
            if (!any(keep)) {
                stop("the mutation table shares no sample identifiers with ",
                     "the expression data", call. = FALSE)
            }
            mut <- mutDF[keep, , drop = FALSE]
        }
    }

    OncoSeqExperiment(
        assays = list(counts = counts),
        colData = cd, rowData = rd,
        mutations = mut, cnv = cnvMat,
        provenance = list(source = "local", call = match.call()))
}

#' Import transcript-level quantifications through tximport
#'
#' @description
#' Summarises transcript-level abundance files produced by salmon, kallisto,
#' RSEM or similar to the gene level with
#' [tximport::tximport()] and returns an \linkS4class{OncoSeqExperiment}
#' holding the gene-level counts, the average transcript lengths and the TPM
#' estimates.
#'
#' @param files Named character vector of quantification files; the names
#'   become sample identifiers.
#' @param tx2gene A two-column `data.frame` mapping transcript to gene
#'   identifiers.
#' @param type Quantification tool, passed to [tximport::tximport()].
#' @param colData Optional sample annotation, matched to `names(files)`.
#' @param countsFromAbundance Passed to [tximport::tximport()]; the default
#'   `"lengthScaledTPM"` yields counts suitable for direct use in
#'   differential expression.
#' @param ... Further arguments passed to [tximport::tximport()].
#'
#' @return An \linkS4class{OncoSeqExperiment} with assays `counts`,
#'   `abundance` and `length`.
#'
#' @examples
#' # Requires the tximport package and quantification files on disk:
#' if (requireNamespace("tximport", quietly = TRUE)) {
#'   # files <- c(S1 = "S1/quant.sf", S2 = "S2/quant.sf")
#'   # ose <- importTximport(files, tx2gene = t2g, type = "salmon")
#' }
#'
#' @export
importTximport <- function(files, tx2gene, type = "salmon", colData = NULL,
                           countsFromAbundance = "lengthScaledTPM", ...) {
    .requireNamespaces("tximport", "importTximport()")
    if (is.null(names(files))) {
        stop("'files' must be a named vector; the names become sample ",
             "identifiers", call. = FALSE)
    }
    missing_files <- files[!file.exists(files)]
    if (length(missing_files) > 0L) {
        stop(length(missing_files), " quantification file(s) do not exist, ",
             "e.g. ", missing_files[1L], call. = FALSE)
    }
    txi <- tximport::tximport(files, type = type, tx2gene = tx2gene,
                              countsFromAbundance = countsFromAbundance, ...)
    assayList <- list(
        counts = round(txi$counts),
        abundance = txi$abundance,
        length = txi$length)
    cd <- if (is.null(colData)) NULL else as.data.frame(colData)
    if (!is.null(cd) && is.null(rownames(cd))) rownames(cd) <- names(files)
    OncoSeqExperiment(assays = assayList, colData = cd,
        provenance = list(source = "tximport", type = type,
                          countsFromAbundance = countsFromAbundance))
}

#' Convert an existing object to an OncoSeqExperiment
#'
#' @description
#' Wraps a \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}, a `DESeqDataSet` or an
#' `edgeR` `DGEList` in an \linkS4class{OncoSeqExperiment} so that it can
#' enter the OncoSeq workflow, optionally attaching mutation and copy-number
#' layers.
#'
#' @param object A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}, `DESeqDataSet` or
#'   `DGEList`.
#' @param mutations Optional variant table.
#' @param cnv Optional gene-level copy-number matrix or
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}.
#'
#' @return An \linkS4class{OncoSeqExperiment}.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = matrix(rpois(100, 30), nrow = 10,
#'         dimnames = list(paste0("G", 1:10), paste0("S", 1:10)))))
#' asOncoSeq(se)
#'
#' @export
asOncoSeq <- function(object, mutations = NULL, cnv = NULL) {
    if (is(object, "OncoSeqExperiment")) {
        if (!is.null(mutations)) mutations(object) <- mutations
        if (!is.null(cnv)) cnv(object) <- cnv
        return(object)
    }
    if (is(object, "DGEList")) {
        cd <- object$samples
        se <- SummarizedExperiment::SummarizedExperiment(
            assays = list(counts = object$counts),
            colData = .asDataFrame(cd))
    } else if (is(object, "SummarizedExperiment")) {
        se <- as(object, "SummarizedExperiment")
    } else {
        stop("'object' must be a SummarizedExperiment, DESeqDataSet, ",
             "DGEList or OncoSeqExperiment; got ", class(object)[1L],
             call. = FALSE)
    }
    OncoSeqExperiment(se, mutations = mutations, cnv = cnv,
        provenance = list(source = paste0("coerced from ", class(object)[1L])))
}

## Read a delimited file when given a path, otherwise return unchanged.
.readIfPath <- function(x, sep = "\t", rowNames = FALSE) {
    if (is.null(x)) return(NULL)
    if (!is.character(x) || length(x) != 1L) return(x)
    if (!file.exists(x)) {
        stop("file not found: ", x, call. = FALSE)
    }
    df <- utils::read.delim(x, sep = sep, header = TRUE, check.names = FALSE,
                            stringsAsFactors = FALSE)
    if (rowNames && ncol(df) > 1L) {
        ids <- as.character(df[[1L]])
        if (anyDuplicated(ids) > 0L) {
            stop("the first column of ", x, " contains duplicated ",
                 "identifiers and cannot be used as row names", call. = FALSE)
        }
        df <- df[, -1L, drop = FALSE]
        rownames(df) <- ids
    }
    df
}
