#' Normalise expression data
#'
#' @description
#' Computes a normalised expression assay and adds it to the object, leaving
#' the raw counts untouched. Count-scale methods wrap the standard
#' implementations in `edgeR` and `DESeq2`.
#'
#' @details
#' Supported methods:
#' \describe{
#'   \item{`"cpm"`}{Counts per million.}
#'   \item{`"logcpm"`}{log2 counts per million with a prior count.}
#'   \item{`"tmm"`}{log2 CPM after trimmed mean of M-values scaling
#'     ([edgeR::normLibSizes()]).}
#'   \item{`"upperquartile"`, `"rle"`}{log2 CPM after the corresponding
#'     `edgeR` scaling factors.}
#'   \item{`"vst"`}{Variance-stabilising transformation
#'     ([DESeq2::vst()]); the standard choice for clustering and
#'     principal-component analysis of count data.}
#'   \item{`"rlog"`}{Regularised log ([DESeq2::rlog()]); better behaved than
#'     `vst` for small experiments but much slower.}
#'   \item{`"log2tpm"`}{log2 transcripts per million, which requires a
#'     feature-length column.}
#' }
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param method Normalisation method; see Details.
#' @param assay Name of the count assay to transform.
#' @param name Name under which the new assay is stored; defaults to the
#'   method name.
#' @param priorCount Prior count for the log transformations.
#' @param lengthColumn Name of a `rowData()` column holding feature lengths
#'   in bases, required for `method = "log2tpm"`.
#' @param blind Passed to [DESeq2::vst()] / [DESeq2::rlog()].
#'
#' @return The \linkS4class{OncoSeqExperiment} with an added assay.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose <- normalizeCounts(ose, method = "logcpm")
#' SummarizedExperiment::assayNames(ose)
#'
#' @export
normalizeCounts <- function(object,
                            method = c("logcpm", "cpm", "tmm",
                                       "upperquartile", "rle", "vst",
                                       "rlog", "log2tpm"),
                            assay = NULL, name = NULL, priorCount = 1,
                            lengthColumn = NULL, blind = TRUE) {
    stopifnot(is(object, "OncoSeqExperiment"))
    method <- match.arg(method)
    counts <- as.matrix(.getAssay(object, assay))
    if (method != "log2tpm") .checkCountsLike(counts, method)

    normalised <- switch(method,
        cpm = edgeR::cpm(counts),
        logcpm = edgeR::cpm(counts, log = TRUE, prior.count = priorCount),
        tmm = .scaledLogCPM(counts, "TMM", priorCount),
        upperquartile = .scaledLogCPM(counts, "upperquartile", priorCount),
        rle = .scaledLogCPM(counts, "RLE", priorCount),
        vst = .deseqTransform(counts, object, "vst", blind),
        rlog = .deseqTransform(counts, object, "rlog", blind),
        log2tpm = .log2TPM(counts, object, lengthColumn, priorCount))

    name <- if (is.null(name)) method else name
    SummarizedExperiment::assay(object, name) <- normalised
    .addProvenance(object, "normalizeCounts", method = method, assay = name)
}

.checkCountsLike <- function(counts, method) {
    if (any(counts < 0, na.rm = TRUE)) {
        stop("method ", sQuote(method), " expects non-negative counts but ",
             "the assay contains negative values. Normalised or ",
             "z-scored data cannot be used here.", call. = FALSE)
    }
    finite <- counts[is.finite(counts)]
    if (length(finite) > 0L &&
        !isTRUE(all.equal(finite, round(finite), tolerance = 1e-8))) {
        warning("method ", sQuote(method), " expects integer counts but the ",
                "assay holds non-integer values; results may be unreliable",
                call. = FALSE)
    }
    invisible(TRUE)
}

.scaledLogCPM <- function(counts, normMethod, priorCount) {
    dge <- edgeR::DGEList(counts = counts)
    dge <- edgeR::normLibSizes(dge, method = normMethod)
    edgeR::cpm(dge, log = TRUE, prior.count = priorCount)
}

.deseqTransform <- function(counts, object, which, blind) {
    storage.mode(counts) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData = SummarizedExperiment::colData(object),
        design = ~ 1)
    res <- if (which == "vst") {
        if (nrow(dds) < 1000L) {
            DESeq2::varianceStabilizingTransformation(dds, blind = blind)
        } else {
            DESeq2::vst(dds, blind = blind)
        }
    } else {
        DESeq2::rlog(dds, blind = blind)
    }
    SummarizedExperiment::assay(res)
}

.log2TPM <- function(counts, object, lengthColumn, priorCount) {
    rd <- SummarizedExperiment::rowData(object)
    if (is.null(lengthColumn)) {
        candidates <- intersect(c("length", "gene_length", "width",
                                  "effective_length"), colnames(rd))
        if (length(candidates) == 0L) {
            stop("method 'log2tpm' needs feature lengths: supply ",
                 "'lengthColumn' naming a rowData column", call. = FALSE)
        }
        lengthColumn <- candidates[1L]
    }
    if (!lengthColumn %in% colnames(rd)) {
        stop("rowData has no column ", sQuote(lengthColumn), call. = FALSE)
    }
    len <- as.numeric(rd[[lengthColumn]])
    if (anyNA(len) || any(len <= 0)) {
        stop("feature lengths must be positive and complete", call. = FALSE)
    }
    rate <- counts / (len / 1000)
    tpm <- t(t(rate) / colSums(rate, na.rm = TRUE)) * 1e6
    log2(tpm + priorCount)
}
