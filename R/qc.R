#' @importFrom matrixStats colMedians rowVars
#' @importFrom stats prcomp mad median quantile cor
NULL

#' Compute per-sample quality-control metrics
#'
#' @description
#' Summarises sequencing depth, detection and expression distribution for
#' every sample and, optionally, the fraction of reads assigned to
#' mitochondrial and ribosomal genes. Metrics are written into `colData()`
#' so that they can be used as covariates or filters downstream.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param assay Name of the assay to summarise; defaults to `counts`.
#' @param detectionThreshold A feature counts as detected in a sample when
#'   its value exceeds this threshold.
#' @param mitoPattern Regular expression identifying mitochondrial genes in
#'   the row names.
#' @param riboPattern Regular expression identifying ribosomal protein genes.
#' @param store When `TRUE` the metrics are added to `colData()` and the
#'   object is returned; when `FALSE` the metric table is returned.
#'
#' @return An \linkS4class{OncoSeqExperiment} with QC columns in `colData()`,
#'   or a \link[S4Vectors:DataFrame-class]{DataFrame} of metrics when `store = FALSE`.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose <- qcMetrics(ose)
#' head(SummarizedExperiment::colData(ose)[, c("library_size",
#'                                             "detected_genes")])
#'
#' @seealso [filterLowExpression()], [detectOutliers()]
#' @export
qcMetrics <- function(object, assay = NULL, detectionThreshold = 0,
                      mitoPattern = "^MT-", riboPattern = "^RP[SL]",
                      store = TRUE) {
    stopifnot(is(object, "OncoSeqExperiment"))
    mat <- .getAssay(object, assay)

    lib <- colSums(mat, na.rm = TRUE)
    detected <- colSums(mat > detectionThreshold, na.rm = TRUE)
    metrics <- DataFrame(
        library_size = lib,
        detected_genes = detected,
        median_expression = matrixStats::colMedians(as.matrix(mat),
                                                    na.rm = TRUE),
        row.names = colnames(object))

    mito <- grepl(mitoPattern, rownames(mat))
    if (any(mito)) {
        metrics$mito_fraction <- colSums(mat[mito, , drop = FALSE],
                                         na.rm = TRUE) / lib
    }
    ribo <- grepl(riboPattern, rownames(mat))
    if (any(ribo)) {
        metrics$ribo_fraction <- colSums(mat[ribo, , drop = FALSE],
                                         na.rm = TRUE) / lib
    }

    if (!store) return(metrics)

    cd <- SummarizedExperiment::colData(object)
    for (nm in colnames(metrics)) cd[[nm]] <- metrics[[nm]]
    SummarizedExperiment::colData(object) <- cd
    object <- .addProvenance(object, "qcMetrics",
                             assay = assay, threshold = detectionThreshold)
    object
}

#' Remove weakly expressed features
#'
#' @description
#' Drops features that carry too little signal to be tested. With a grouping
#' variable the filter is [edgeR::filterByExpr()], the standard choice for
#' count data, which keeps features with enough counts in a number of
#' samples equal to the smallest group. Without a grouping variable a simple
#' counts-per-million rule is applied.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param group Name of a `colData()` column giving the experimental groups,
#'   or a vector of the same length as the number of samples. `NULL` applies
#'   the CPM rule.
#' @param assay Name of the count assay; defaults to `counts`.
#' @param minCount Minimum count required, passed to
#'   [edgeR::filterByExpr()].
#' @param minCPM CPM threshold used when `group` is `NULL`.
#' @param minSamples Number of samples that must exceed `minCPM` when
#'   `group` is `NULL`; defaults to the size of the smallest conceivable
#'   group, 2.
#' @param verbose Report how many features were removed.
#'
#' @return The filtered \linkS4class{OncoSeqExperiment}.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' dim(ose)
#' dim(filterLowExpression(ose, group = "condition"))
#'
#' @export
filterLowExpression <- function(object, group = NULL, assay = NULL,
                                minCount = 10, minCPM = 1, minSamples = 2L,
                                verbose = TRUE) {
    stopifnot(is(object, "OncoSeqExperiment"))
    mat <- as.matrix(.getAssay(object, assay))

    if (!is.null(group)) {
        grp <- .resolveColumn(object, group, "group")
        keep <- edgeR::filterByExpr(mat, group = grp, min.count = minCount)
    } else {
        cpm <- edgeR::cpm(mat)
        keep <- rowSums(cpm > minCPM) >= minSamples
    }
    if (verbose) {
        message("filterLowExpression: keeping ", sum(keep), " of ",
                length(keep), " features")
    }
    out <- object[keep, ]
    .addProvenance(out, "filterLowExpression", group = group,
                   minCount = minCount, kept = sum(keep),
                   removed = sum(!keep))
}

#' Flag candidate outlier samples
#'
#' @description
#' Identifies samples that sit far from the bulk of the data in two
#' complementary ways: distance from the centroid in the space of the
#' leading principal components, and low median correlation with all other
#' samples. Both criteria are reported; nothing is removed.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param assay Assay to use. When absent, log2 counts-per-million are
#'   computed from the counts.
#' @param nTopGenes Number of most variable features used for the principal
#'   component analysis.
#' @param nPC Number of principal components used for the distance.
#' @param madThreshold A sample is flagged when its centroid distance exceeds
#'   the median distance by this many median absolute deviations.
#' @param corMethod Correlation method for the second criterion.
#'
#' @return A \link[S4Vectors:DataFrame-class]{DataFrame} with one row per sample giving the
#'   centroid distance, the median inter-sample correlation and a logical
#'   `outlier` column.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' out <- detectOutliers(ose)
#' sum(out$outlier)
#'
#' @export
detectOutliers <- function(object, assay = NULL, nTopGenes = 500L, nPC = 2L,
                           madThreshold = 3, corMethod = "spearman") {
    stopifnot(is(object, "OncoSeqExperiment"))
    mat <- .logExpressionMatrix(object, assay)
    if (ncol(mat) < 3L) {
        stop("outlier detection needs at least three samples", call. = FALSE)
    }

    vars <- matrixStats::rowVars(mat, na.rm = TRUE)
    top <- head(order(vars, decreasing = TRUE),
                min(nTopGenes, nrow(mat)))
    sub <- mat[top, , drop = FALSE]

    pca <- stats::prcomp(t(sub), scale. = FALSE)
    nPC <- min(nPC, ncol(pca$x))
    coords <- pca$x[, seq_len(nPC), drop = FALSE]
    centre <- colMeans(coords)
    dist <- sqrt(rowSums(sweep(coords, 2, centre, "-")^2))

    cmat <- stats::cor(sub, method = corMethod, use = "pairwise.complete.obs")
    diag(cmat) <- NA_real_
    medCor <- apply(cmat, 2, stats::median, na.rm = TRUE)

    distCut <- stats::median(dist) + madThreshold * stats::mad(dist)
    corCut <- stats::median(medCor) - madThreshold * stats::mad(medCor)

    DataFrame(
        pc_distance = dist,
        median_correlation = medCor,
        distance_outlier = dist > distCut,
        correlation_outlier = medCor < corCut,
        outlier = dist > distCut | medCor < corCut,
        row.names = colnames(object))
}

## Resolve a colData column name or a supplied vector to a factor/vector.
.resolveColumn <- function(object, x, argName) {
    if (length(x) == 1L && is.character(x)) {
        cd <- SummarizedExperiment::colData(object)
        if (!x %in% colnames(cd)) {
            stop(argName, " = ", sQuote(x),
                 " is not a column of colData(object). Available: ",
                 paste(colnames(cd), collapse = ", "), call. = FALSE)
        }
        return(cd[[x]])
    }
    if (length(x) != ncol(object)) {
        stop(argName, " must name a colData column or be a vector of length ",
             ncol(object), call. = FALSE)
    }
    x
}

## Return a log-scale expression matrix, deriving log2-CPM from counts when
## no transformed assay is available.
.logExpressionMatrix <- function(object, assay = NULL) {
    nms <- SummarizedExperiment::assayNames(object)
    if (!is.null(assay)) return(as.matrix(.getAssay(object, assay)))
    for (candidate in c("vst", "rlog", "logcpm", "log2tpm")) {
        if (candidate %in% nms) return(as.matrix(.getAssay(object, candidate)))
    }
    if ("counts" %in% nms) {
        return(edgeR::cpm(as.matrix(.getAssay(object, "counts")),
                          log = TRUE, prior.count = 1))
    }
    mat <- as.matrix(.getAssay(object, nms[1L]))
    if (max(mat, na.rm = TRUE) > 100) log2(mat + 1) else mat
}
