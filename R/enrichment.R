#' @importFrom fgsea fgseaMultilevel
#' @importFrom BiocParallel SerialParam
#' @importFrom stats phyper
NULL

#' Gene set enrichment analysis on a differential-expression result
#'
#' @description
#' Runs [fgsea::fgseaMultilevel()] on features ranked by a statistic taken
#' from a stored differential-expression table. The multilevel algorithm is
#' the current fgsea default and gives accurate p-values far into the tail.
#'
#' @param object An \linkS4class{OncoSeqExperiment} holding a
#'   differential-expression result.
#' @param geneSets A named list of gene identifiers; see [readGMT()].
#' @param deName Name of the stored differential-expression result.
#' @param rankBy Statistic used to rank features. `"signed_log_p"` is
#'   `sign(log2FC) * -log10(pvalue)`, which separates direction from
#'   magnitude better than the fold change alone; `"log2FC"` and `"stat"`
#'   are also available.
#' @param minSize,maxSize Gene-set size limits.
#' @param name Name under which the result is stored.
#' @param BPPARAM A \link[BiocParallel:BiocParallelParam-class]{BiocParallelParam} instance controlling
#'   parallelisation. The default is serial evaluation, which is
#'   reproducible and avoids the socket set-up that fails in sandboxed
#'   environments; pass [BiocParallel::MulticoreParam()] to use several
#'   cores.
#' @param ... Further arguments passed to [fgsea::fgseaMultilevel()].
#'
#' @return The object with a result table stored in the registry, carrying
#'   columns `set`, `size`, `ES`, `NES`, `pvalue`, `padj`, `leadingEdge`
#'   and `method`, ordered by adjusted p-value.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose <- diffExpress(ose, ~ condition,
#'                    contrast = c("condition", "tumour", "normal"),
#'                    engine = "edgeR", verbose = FALSE)
#' sets <- readGMT(system.file("extdata", "example_sets.gmt",
#'                             package = "OncoSeq"))
#' ose <- runGSEA(ose, sets)
#' head(getAnalysis(ose, "gsea")[, c("set", "NES", "padj")])
#'
#' @seealso [runORA()], [plotEnrichment()]
#' @export
runGSEA <- function(object, geneSets, deName = "de",
                    rankBy = c("signed_log_p", "stat", "log2FC"),
                    minSize = 5L, maxSize = 500L, name = "gsea",
                    BPPARAM = BiocParallel::SerialParam(), ...) {
    stopifnot(is(object, "OncoSeqExperiment"))
    rankBy <- match.arg(rankBy)
    de <- getAnalysis(object, deName)
    geneSets <- .checkGeneSets(geneSets, de$feature)

    ranks <- .rankStatistic(de, rankBy)
    res <- fgsea::fgseaMultilevel(pathways = geneSets, stats = ranks,
                                  minSize = minSize, maxSize = maxSize,
                                  BPPARAM = BPPARAM, ...)
    res <- as.data.frame(res)
    out <- data.frame(
        set = res$pathway,
        size = res$size,
        ES = res$ES,
        NES = res$NES,
        pvalue = res$pval,
        padj = res$padj,
        stringsAsFactors = FALSE)
    out$leadingEdge <- vapply(res$leadingEdge, paste, character(1),
                              collapse = ";")
    out$method <- "fgsea"
    out <- out[order(out$padj, out$pvalue, na.last = TRUE), , drop = FALSE]
    rownames(out) <- NULL

    object <- .addAnalysis(object, name, out, analysis = "gene set enrichment",
        deName = deName, rankBy = rankBy, nSets = length(geneSets))
    .addProvenance(object, "runGSEA", name = name, rankBy = rankBy)
}

#' Over-representation analysis of a differentially expressed gene list
#'
#' @description
#' Tests whether each gene set is over-represented among the significant
#' features of a stored differential-expression result, using the
#' hypergeometric distribution with the tested features as the universe and
#' Benjamini-Hochberg correction across sets.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param geneSets A named list of gene identifiers.
#' @param deName Name of the stored differential-expression result.
#' @param alpha Adjusted-p-value cut-off defining the significant list.
#' @param lfcThreshold Minimum absolute log2 fold change for inclusion.
#' @param direction Restrict the list to `"up"`, `"down"` or `"both"`.
#' @param universe Character vector of background features; defaults to all
#'   features in the differential-expression table, which is the correct
#'   universe because only those were testable.
#' @param minSize,maxSize Gene-set size limits, applied after intersecting
#'   each set with the universe.
#' @param name Name under which the result is stored.
#'
#' @return The object with a result table in the registry, carrying `set`,
#'   `size`, `overlap`, `expected`, `foldEnrichment`, `pvalue`, `padj`,
#'   `genes` and `method`.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose <- diffExpress(ose, ~ condition,
#'                    contrast = c("condition", "tumour", "normal"),
#'                    engine = "edgeR", verbose = FALSE)
#' sets <- readGMT(system.file("extdata", "example_sets.gmt",
#'                             package = "OncoSeq"))
#' ose <- runORA(ose, sets, alpha = 0.1)
#' head(getAnalysis(ose, "ora")[, c("set", "overlap", "padj")])
#'
#' @seealso [runGSEA()]
#' @export
runORA <- function(object, geneSets, deName = "de", alpha = 0.05,
                   lfcThreshold = 0, direction = c("both", "up", "down"),
                   universe = NULL, minSize = 5L, maxSize = 500L,
                   name = "ora") {
    stopifnot(is(object, "OncoSeqExperiment"))
    direction <- match.arg(direction)
    de <- getAnalysis(object, deName)

    if (is.null(universe)) universe <- de$feature
    universe <- unique(as.character(universe))

    sig <- de$feature[which(!is.na(de$padj) & de$padj < alpha &
                            abs(de$log2FC) >= lfcThreshold)]
    if (direction == "up") {
        sig <- de$feature[which(!is.na(de$padj) & de$padj < alpha &
                                de$log2FC >= lfcThreshold)]
    } else if (direction == "down") {
        sig <- de$feature[which(!is.na(de$padj) & de$padj < alpha &
                                de$log2FC <= -lfcThreshold)]
    }
    sig <- intersect(unique(sig), universe)
    if (length(sig) == 0L) {
        stop("no features pass alpha = ", alpha, " and |log2FC| >= ",
             lfcThreshold, "; relax the thresholds", call. = FALSE)
    }

    geneSets <- .checkGeneSets(geneSets, universe)
    sets <- lapply(geneSets, intersect, universe)
    sizes <- lengths(sets)
    sets <- sets[sizes >= minSize & sizes <= maxSize]
    if (length(sets) == 0L) {
        stop("no gene set has between ", minSize, " and ", maxSize,
             " members in the universe", call. = FALSE)
    }

    N <- length(universe)
    k <- length(sig)
    overlapGenes <- lapply(sets, intersect, sig)
    q <- lengths(overlapGenes)
    m <- lengths(sets)
    ## P(X >= q) with X ~ Hypergeometric(m white, N - m black, k drawn)
    pvalue <- stats::phyper(q - 1L, m, N - m, k, lower.tail = FALSE)
    expected <- k * m / N

    out <- data.frame(
        set = names(sets),
        size = as.integer(m),
        overlap = as.integer(q),
        expected = expected,
        foldEnrichment = ifelse(expected > 0, q / expected, NA_real_),
        pvalue = pvalue,
        padj = .bh(pvalue),
        genes = vapply(overlapGenes, paste, character(1), collapse = ";"),
        method = "hypergeometric",
        stringsAsFactors = FALSE)
    out <- out[order(out$padj, out$pvalue, na.last = TRUE), , drop = FALSE]
    rownames(out) <- NULL

    object <- .addAnalysis(object, name, out,
        analysis = "over-representation", deName = deName, alpha = alpha,
        direction = direction, nSignificant = k, universe = N)
    .addProvenance(object, "runORA", name = name, direction = direction)
}

.checkGeneSets <- function(geneSets, features) {
    if (!is.list(geneSets) || is.null(names(geneSets))) {
        stop("'geneSets' must be a named list of character vectors; see ",
             "readGMT()", call. = FALSE)
    }
    overlap <- length(intersect(unique(unlist(geneSets, use.names = FALSE)),
                                features))
    if (overlap == 0L) {
        stop("no gene-set member matches any feature identifier. The gene ",
             "sets use ", paste(utils::head(unlist(geneSets,
                use.names = FALSE), 3L), collapse = ", "),
             " while the data use ", paste(utils::head(features, 3L),
                collapse = ", "),
             ". Map identifiers to a common namespace first.",
             call. = FALSE)
    }
    geneSets
}

.rankStatistic <- function(de, rankBy) {
    stat <- switch(rankBy,
        log2FC = de$log2FC,
        stat = de$stat,
        signed_log_p = sign(de$log2FC) * -log10(de$pvalue))
    if (all(is.na(stat))) {
        stop("ranking statistic ", sQuote(rankBy), " is entirely missing in ",
             "this result table; try rankBy = 'log2FC'", call. = FALSE)
    }
    names(stat) <- de$feature
    stat <- stat[is.finite(stat)]
    stat <- stat[!duplicated(names(stat))]
    sort(stat, decreasing = TRUE)
}
