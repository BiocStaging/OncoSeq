#' Summarise gene-level copy-number alterations
#'
#' @description
#' Counts amplifications and deletions per gene and per sample from the
#' discrete copy-number layer. Thresholds follow the GISTIC convention in
#' which `2` and `-2` denote high-level amplification and deep deletion and
#' `1`/`-1` denote single-copy gain and loss.
#'
#' @param object An \linkS4class{OncoSeqExperiment} with a copy-number layer.
#' @param ampThreshold Values greater than or equal to this count as a gain.
#' @param delThreshold Values less than or equal to this count as a loss.
#' @param assay Name or index of the copy-number assay.
#'
#' @return A list with `perGene` and `perSample` `data.frame`s.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' cs <- cnvSummary(ose)
#' head(cs$perGene)
#'
#' @export
cnvSummary <- function(object, ampThreshold = 2, delThreshold = -2,
                       assay = 1L) {
    mat <- .cnvMatrix(object, assay)
    gain <- mat >= ampThreshold
    loss <- mat <= delThreshold

    perGene <- data.frame(
        gene = rownames(mat),
        n_gain = rowSums(gain, na.rm = TRUE),
        n_loss = rowSums(loss, na.rm = TRUE),
        stringsAsFactors = FALSE)
    perGene$n_samples <- ncol(mat)
    perGene$freq_gain <- perGene$n_gain / ncol(mat)
    perGene$freq_loss <- perGene$n_loss / ncol(mat)
    perGene <- perGene[order(-pmax(perGene$freq_gain, perGene$freq_loss)), ,
                       drop = FALSE]
    rownames(perGene) <- NULL

    perSample <- data.frame(
        sample = colnames(mat),
        n_gain = colSums(gain, na.rm = TRUE),
        n_loss = colSums(loss, na.rm = TRUE),
        stringsAsFactors = FALSE)
    perSample$fraction_altered <-
        (perSample$n_gain + perSample$n_loss) / nrow(mat)
    rownames(perSample) <- NULL

    list(perGene = perGene, perSample = perSample)
}

#' Copy-number alteration frequency
#'
#' @description
#' Per-gene frequency of gain and loss across the cohort, in the long format
#' expected by [plotCNVFrequency()].
#'
#' @inheritParams cnvSummary
#' @param genes Optional character vector restricting the output.
#' @param topN Keep the `topN` genes with the highest alteration frequency.
#'
#' @return A `data.frame` with columns `gene`, `alteration` (`"gain"` or
#'   `"loss"`) and `frequency`.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' head(cnvFrequency(ose, topN = 5))
#'
#' @export
cnvFrequency <- function(object, genes = NULL, topN = 20, ampThreshold = 2,
                         delThreshold = -2, assay = 1L) {
    perGene <- cnvSummary(object, ampThreshold, delThreshold, assay)$perGene
    if (!is.null(genes)) {
        perGene <- perGene[perGene$gene %in% genes, , drop = FALSE]
    }
    if (is.finite(topN)) perGene <- utils::head(perGene, topN)
    rbind(
        data.frame(gene = perGene$gene, alteration = "gain",
                   frequency = perGene$freq_gain, stringsAsFactors = FALSE),
        data.frame(gene = perGene$gene, alteration = "loss",
                   frequency = perGene$freq_loss, stringsAsFactors = FALSE))
}

#' Correlate copy number with expression
#'
#' @description
#' Tests the dosage effect of copy-number alteration on transcript abundance
#' gene by gene, over the samples that carry both layers. A positive
#' correlation is the expected signature of a copy-number-driven expression
#' change and is used to separate passenger from dosage-sensitive genes.
#'
#' @param object An \linkS4class{OncoSeqExperiment} with expression and
#'   copy-number layers.
#' @param genes Genes to test; defaults to every gene present in both
#'   layers.
#' @param exprAssay Expression assay to use; log-scale values are preferred
#'   and are derived from counts when necessary.
#' @param cnvAssay Name or index of the copy-number assay.
#' @param method Correlation method, passed to [stats::cor.test()].
#' @param minSamples Minimum number of complete pairs required to test a
#'   gene.
#'
#' @return A `data.frame` with `gene`, `n`, `estimate`, `pvalue` and `padj`,
#'   ordered by adjusted p-value.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' head(cnvExpressionCorrelation(ose))
#'
#' @export
cnvExpressionCorrelation <- function(object, genes = NULL, exprAssay = NULL,
                                     cnvAssay = 1L, method = "spearman",
                                     minSamples = 10L) {
    cnvMat <- .cnvMatrix(object, cnvAssay)
    exprMat <- .logExpressionMatrix(object, exprAssay)

    commonSamples <- intersect(colnames(exprMat), colnames(cnvMat))
    if (length(commonSamples) < minSamples) {
        stop("only ", length(commonSamples), " sample(s) carry both ",
             "expression and copy-number data; at least ", minSamples,
             " are needed", call. = FALSE)
    }
    commonGenes <- intersect(rownames(exprMat), rownames(cnvMat))
    if (!is.null(genes)) commonGenes <- intersect(commonGenes, genes)
    if (length(commonGenes) == 0L) {
        stop("no gene identifier is shared by the expression and ",
             "copy-number layers", call. = FALSE)
    }

    e <- exprMat[commonGenes, commonSamples, drop = FALSE]
    c0 <- cnvMat[commonGenes, commonSamples, drop = FALSE]

    res <- lapply(seq_along(commonGenes), function(i) {
        x <- as.numeric(c0[i, ])
        y <- as.numeric(e[i, ])
        ok <- is.finite(x) & is.finite(y)
        if (sum(ok) < minSamples || length(unique(x[ok])) < 2L) {
            return(c(n = sum(ok), estimate = NA_real_, pvalue = NA_real_))
        }
        tt <- suppressWarnings(
            stats::cor.test(x[ok], y[ok], method = method))
        c(n = sum(ok), estimate = unname(tt$estimate),
          pvalue = unname(tt$p.value))
    })
    out <- as.data.frame(do.call(rbind, res))
    out$gene <- commonGenes
    out$padj <- .bh(out$pvalue)
    out <- out[order(out$padj, out$pvalue, na.last = TRUE),
               c("gene", "n", "estimate", "pvalue", "padj"), drop = FALSE]
    rownames(out) <- NULL
    out
}

.cnvMatrix <- function(object, assay = 1L) {
    stopifnot(is(object, "OncoSeqExperiment"))
    cn <- cnv(object)
    if (is.null(cn)) {
        stop("this object carries no copy-number data; attach it with ",
             "cnv(object) <- ...", call. = FALSE)
    }
    as.matrix(SummarizedExperiment::assay(cn, assay))
}
