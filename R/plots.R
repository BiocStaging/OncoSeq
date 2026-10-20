#' @importFrom ggplot2 ggplot aes geom_point geom_text geom_tile geom_col
#'   geom_step geom_hline geom_vline geom_errorbar labs theme_bw theme
#'   element_text element_blank scale_colour_manual scale_fill_manual
#'   scale_fill_gradient2 coord_flip scale_y_continuous .data
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices colorRampPalette
NULL

.oncoTheme <- function(base_size = 11) {
    ggplot2::theme_bw(base_size = base_size) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            strip.background = ggplot2::element_blank(),
            legend.key = ggplot2::element_blank())
}

#' Principal component plot of samples
#'
#' @description
#' Projects samples onto their leading principal components, computed from
#' the most variable features of a log-scale assay, and colours them by a
#' sample annotation.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param colourBy `colData()` column used for point colour.
#' @param shapeBy Optional `colData()` column used for point shape.
#' @param assay Assay to use; log2 counts-per-million are derived from
#'   counts when no transformed assay is present.
#' @param nTopGenes Number of most variable features used.
#' @param components Length-two integer vector giving which components to
#'   draw.
#' @param label Draw sample names next to the points.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' plotSamplePCA(ose, colourBy = "condition")
#'
#' @export
plotSamplePCA <- function(object, colourBy = NULL, shapeBy = NULL,
                          assay = NULL, nTopGenes = 500L,
                          components = c(1L, 2L), label = FALSE) {
    mat <- .logExpressionMatrix(object, assay)
    vars <- matrixStats::rowVars(mat, na.rm = TRUE)
    top <- utils::head(order(vars, decreasing = TRUE),
                       min(nTopGenes, nrow(mat)))
    pca <- stats::prcomp(t(mat[top, , drop = FALSE]), scale. = FALSE)
    pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

    df <- data.frame(
        PCx = pca$x[, components[1L]],
        PCy = pca$x[, components[2L]],
        sample = colnames(object), stringsAsFactors = FALSE)
    cd <- as.data.frame(SummarizedExperiment::colData(object))
    if (!is.null(colourBy)) df$colour <- .resolveColumn(object, colourBy,
                                                        "colourBy")
    if (!is.null(shapeBy)) df$shape <- factor(.resolveColumn(object, shapeBy,
                                                             "shapeBy"))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$PCx, y = .data$PCy))
    pointMapping <- if (!is.null(colourBy) && !is.null(shapeBy)) {
        ggplot2::aes(colour = .data$colour, shape = .data$shape)
    } else if (!is.null(colourBy)) {
        ggplot2::aes(colour = .data$colour)
    } else if (!is.null(shapeBy)) {
        ggplot2::aes(shape = .data$shape)
    } else {
        ggplot2::aes()
    }
    p <- p + ggplot2::geom_point(mapping = pointMapping, size = 2.4,
                                 alpha = 0.9)
    if (label) {
        p <- p + ggrepel::geom_text_repel(
            ggplot2::aes(label = .data$sample), size = 2.8,
            max.overlaps = 20)
    }
    p + ggplot2::labs(
            x = paste0("PC", components[1L], " (", pct[components[1L]], "%)"),
            y = paste0("PC", components[2L], " (", pct[components[2L]], "%)"),
            colour = colourBy, shape = shapeBy) +
        .oncoTheme()
}

#' Volcano plot of a differential-expression result
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param deName Name of the stored differential-expression result.
#' @param alpha Adjusted-p-value threshold for calling a feature
#'   significant.
#' @param lfcThreshold Absolute log2 fold-change threshold.
#' @param nLabel Number of top features to label, by adjusted p-value.
#' @param labelGenes Explicit features to label, overriding `nLabel`.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose <- diffExpress(ose, ~ condition,
#'                    contrast = c("condition", "tumour", "normal"),
#'                    engine = "edgeR", verbose = FALSE)
#' plotVolcano(ose)
#'
#' @export
plotVolcano <- function(object, deName = "de", alpha = 0.05,
                        lfcThreshold = 1, nLabel = 10L,
                        labelGenes = NULL) {
    de <- getAnalysis(object, deName)
    de <- de[is.finite(de$pvalue), , drop = FALSE]
    de$negLogP <- -log10(de$pvalue)
    de$status <- "not significant"
    sig <- !is.na(de$padj) & de$padj < alpha & abs(de$log2FC) >= lfcThreshold
    de$status[sig & de$log2FC > 0] <- "up"
    de$status[sig & de$log2FC < 0] <- "down"
    de$status <- factor(de$status, levels = c("down", "not significant", "up"))

    lab <- if (!is.null(labelGenes)) {
        de[de$feature %in% labelGenes, , drop = FALSE]
    } else {
        utils::head(de[sig, , drop = FALSE], nLabel)
    }

    ggplot2::ggplot(de, ggplot2::aes(x = .data$log2FC, y = .data$negLogP)) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$status),
                            size = 1.2, alpha = 0.7) +
        ggplot2::geom_vline(xintercept = c(-lfcThreshold, lfcThreshold),
                            linetype = "dashed", linewidth = 0.3) +
        ggplot2::scale_colour_manual(values = c(
            down = "#2166AC", `not significant` = "grey72", up = "#B2182B")) +
        ggrepel::geom_text_repel(data = lab,
            ggplot2::aes(label = .data$feature), size = 2.9,
            max.overlaps = 30, min.segment.length = 0) +
        ggplot2::labs(x = expression(log[2]~fold~change),
                      y = expression(-log[10]~italic(P)),
                      colour = NULL) +
        .oncoTheme()
}

#' Mean-difference (MA) plot of a differential-expression result
#'
#' @inheritParams plotVolcano
#'
#' @return A `ggplot` object.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose <- diffExpress(ose, ~ condition,
#'                    contrast = c("condition", "tumour", "normal"),
#'                    engine = "edgeR", verbose = FALSE)
#' plotMAplot(ose)
#'
#' @export
plotMAplot <- function(object, deName = "de", alpha = 0.05,
                       lfcThreshold = 0) {
    de <- getAnalysis(object, deName)
    de <- de[is.finite(de$baseMean) & de$baseMean > 0, , drop = FALSE]
    de$significant <- !is.na(de$padj) & de$padj < alpha &
        abs(de$log2FC) >= lfcThreshold

    ggplot2::ggplot(de, ggplot2::aes(x = log2(.data$baseMean),
                                     y = .data$log2FC)) +
        ggplot2::geom_point(ggplot2::aes(colour = .data$significant),
                            size = 1.1, alpha = 0.6) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
        ggplot2::scale_colour_manual(
            values = c(`FALSE` = "grey72", `TRUE` = "#B2182B"),
            labels = c(`FALSE` = "not significant", `TRUE` = "significant")) +
        ggplot2::labs(x = expression(log[2]~mean~expression),
                      y = expression(log[2]~fold~change), colour = NULL) +
        .oncoTheme()
}

#' Heatmap of selected features
#'
#' @description
#' Draws a per-feature standardised (z-scored) heatmap of the most variable
#' or explicitly chosen features, with optional sample annotation shown as a
#' strip above the matrix.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param features Features to draw; defaults to the `nTopGenes` most
#'   variable.
#' @param nTopGenes Number of features when `features` is `NULL`.
#' @param assay Assay to use.
#' @param annotationColumn Optional `colData()` column used to order and
#'   annotate samples.
#' @param scale Standardise each feature across samples.
#' @param clusterSamples Order samples by hierarchical clustering within
#'   annotation groups.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' plotExpressionHeatmap(ose, nTopGenes = 20, annotationColumn = "condition")
#'
#' @export
plotExpressionHeatmap <- function(object, features = NULL, nTopGenes = 50L,
                                  assay = NULL, annotationColumn = NULL,
                                  scale = TRUE, clusterSamples = TRUE) {
    mat <- .logExpressionMatrix(object, assay)
    if (is.null(features)) {
        vars <- matrixStats::rowVars(mat, na.rm = TRUE)
        features <- rownames(mat)[utils::head(order(vars, decreasing = TRUE),
                                              min(nTopGenes, nrow(mat)))]
    }
    missing_features <- setdiff(features, rownames(mat))
    if (length(missing_features) > 0L) {
        warning(length(missing_features), " requested feature(s) are not in ",
                "the assay and were dropped", call. = FALSE)
        features <- intersect(features, rownames(mat))
    }
    sub <- mat[features, , drop = FALSE]
    if (scale) {
        sub <- t(scale(t(sub)))
        sub[!is.finite(sub)] <- 0
    }

    sampleOrder <- colnames(sub)
    if (!is.null(annotationColumn)) {
        ann <- .resolveColumn(object, annotationColumn, "annotationColumn")
        sampleOrder <- colnames(sub)[order(ann)]
    }
    if (clusterSamples && ncol(sub) > 2L) {
        if (is.null(annotationColumn)) {
            hc <- stats::hclust(stats::dist(t(sub)))
            sampleOrder <- colnames(sub)[hc$order]
        } else {
            ann <- .resolveColumn(object, annotationColumn, "annotationColumn")
            sampleOrder <- unlist(lapply(split(colnames(sub), ann),
                function(ss) {
                    if (length(ss) < 3L) return(ss)
                    ss[stats::hclust(stats::dist(t(sub[, ss,
                        drop = FALSE])))$order]
                }), use.names = FALSE)
        }
    }
    featureOrder <- if (nrow(sub) > 2L) {
        rownames(sub)[stats::hclust(stats::dist(sub))$order]
    } else rownames(sub)

    long <- data.frame(
        feature = factor(rep(rownames(sub), times = ncol(sub)),
                         levels = featureOrder),
        sample = factor(rep(colnames(sub), each = nrow(sub)),
                        levels = sampleOrder),
        value = as.vector(sub), stringsAsFactors = FALSE)

    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$sample,
                                            y = .data$feature,
                                            fill = .data$value)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white",
                                      high = "#B2182B", midpoint = 0) +
        ggplot2::labs(x = NULL, y = NULL,
                      fill = if (scale) "z-score" else "expression") +
        .oncoTheme() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_text(size = 6))
    if (!is.null(annotationColumn)) {
        p <- p + ggplot2::labs(subtitle = paste("samples ordered by",
                                                annotationColumn))
    }
    p
}

#' Kaplan-Meier curves
#'
#' @description
#' Draws the survival curves from a [survivalFit()] result, annotated with
#' the log-rank p-value and the number at risk in each group.
#'
#' @param fit An object returned by [survivalFit()].
#' @param showPValue Annotate the panel with the log-rank p-value.
#' @param palette Optional vector of colours, one per group.
#' @param timeLabel Axis label for follow-up time.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' plotKM(survivalFit(ose, "TP53", by = "mutation"))
#'
#' @export
plotKM <- function(fit, showPValue = TRUE, palette = NULL,
                   timeLabel = "Time") {
    if (!inherits(fit, "oncoseq_survival")) {
        stop("'fit' must be the result of survivalFit()", call. = FALSE)
    }
    sf <- fit$fit
    strata <- if (is.null(sf$strata)) {
        rep("all", length(sf$time))
    } else {
        rep(sub("^group=", "", names(sf$strata)), sf$strata)
    }
    df <- data.frame(time = sf$time, surv = sf$surv, group = strata,
                     stringsAsFactors = FALSE)
    starts <- do.call(rbind, lapply(unique(df$group), function(g) {
        data.frame(time = 0, surv = 1, group = g, stringsAsFactors = FALSE)
    }))
    df <- rbind(starts, df)
    df$group <- factor(df$group, levels = fit$groups)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$time, y = .data$surv,
                                          colour = .data$group)) +
        ggplot2::geom_step(linewidth = 0.7) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::labs(x = timeLabel, y = "Survival probability",
                      colour = NULL,
                      title = paste(fit$feature, "-", fit$layer)) +
        .oncoTheme()
    if (!is.null(palette)) {
        p <- p + ggplot2::scale_colour_manual(values = palette)
    }
    if (showPValue) {
        lab <- paste0("log-rank P = ",
                      format.pval(fit$p_logrank, digits = 2, eps = 1e-16))
        p <- p + ggplot2::labs(subtitle = paste0(lab, "   (n = ",
            paste(fit$n, collapse = ", "), ")"))
    }
    p
}

#' Oncoprint-style mutation landscape
#'
#' @description
#' Tile plot of the matrix returned by [oncoplotData()]: genes in rows,
#' samples in columns, coloured by variant classification, with the
#' per-gene mutation frequency on the right.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param genes Genes to show; defaults to the `topN` most frequently
#'   mutated.
#' @param topN Number of genes when `genes` is `NULL`.
#' @param nonsynonymousOnly Restrict to protein-altering variants.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' plotOncoprint(ose, topN = 8)
#'
#' @export
plotOncoprint <- function(object, genes = NULL, topN = 20L,
                          nonsynonymousOnly = TRUE) {
    mat <- oncoplotData(object, genes = genes, topN = topN,
                        nonsynonymousOnly = nonsynonymousOnly)
    freq <- rowSums(mat != "") / ncol(mat)
    geneOrder <- names(sort(freq))

    long <- data.frame(
        gene = factor(rep(rownames(mat), times = ncol(mat)),
                      levels = geneOrder),
        sample = factor(rep(colnames(mat), each = nrow(mat)),
                        levels = colnames(mat)),
        class = as.vector(mat), stringsAsFactors = FALSE)
    long$class[long$class == ""] <- NA_character_

    classes <- sort(unique(stats::na.omit(long$class)))
    cols <- stats::setNames(
        grDevices::colorRampPalette(
            c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
              "#E6AB02"))(length(classes)), classes)

    ggplot2::ggplot(long, ggplot2::aes(x = .data$sample, y = .data$gene)) +
        ggplot2::geom_tile(fill = "grey93", colour = NA) +
        ggplot2::geom_tile(data = long[!is.na(long$class), , drop = FALSE],
                           ggplot2::aes(fill = .data$class)) +
        ggplot2::scale_fill_manual(values = cols, na.value = "grey93") +
        ggplot2::labs(x = paste0(ncol(mat), " samples"), y = NULL,
                      fill = "variant class") +
        .oncoTheme() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                       axis.ticks.x = ggplot2::element_blank())
}

#' Bar chart of enrichment results
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param name Name of a stored [runGSEA()] or [runORA()] result.
#' @param topN Number of gene sets to show.
#' @param alpha Only sets with an adjusted p-value below this are shown;
#'   `NULL` shows the top sets regardless.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose <- diffExpress(ose, ~ condition,
#'                    contrast = c("condition", "tumour", "normal"),
#'                    engine = "edgeR", verbose = FALSE)
#' sets <- readGMT(system.file("extdata", "example_sets.gmt",
#'                             package = "OncoSeq"))
#' ose <- runGSEA(ose, sets)
#' plotEnrichmentBar(ose, "gsea")
#'
#' @export
plotEnrichmentBar <- function(object, name = "gsea", topN = 15L,
                              alpha = NULL) {
    res <- getAnalysis(object, name)
    if (!is.null(alpha)) {
        res <- res[!is.na(res$padj) & res$padj < alpha, , drop = FALSE]
        if (nrow(res) == 0L) {
            stop("no gene set reaches padj < ", alpha, call. = FALSE)
        }
    }
    res <- utils::head(res, topN)
    valueCol <- if ("NES" %in% colnames(res)) "NES" else "foldEnrichment"
    res$set <- factor(res$set, levels = res$set[order(res[[valueCol]])])

    ggplot2::ggplot(res, ggplot2::aes(x = .data$set, y = .data[[valueCol]],
                                      fill = -log10(.data$padj))) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::labs(x = NULL, y = valueCol,
                      fill = expression(-log[10]~adjusted~italic(P))) +
        .oncoTheme() +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7))
}

#' Copy-number alteration frequency plot
#'
#' @param object An \linkS4class{OncoSeqExperiment} with a copy-number layer.
#' @param genes Optional genes to show.
#' @param topN Number of genes with the highest alteration frequency.
#' @param ampThreshold,delThreshold Discrete thresholds for gain and loss.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' plotCNVFrequency(ose, topN = 10)
#'
#' @export
plotCNVFrequency <- function(object, genes = NULL, topN = 20L,
                             ampThreshold = 2, delThreshold = -2) {
    df <- cnvFrequency(object, genes = genes, topN = topN,
                       ampThreshold = ampThreshold,
                       delThreshold = delThreshold)
    df$signed <- ifelse(df$alteration == "loss", -df$frequency, df$frequency)
    order_genes <- unique(df$gene[order(df$signed)])
    df$gene <- factor(df$gene, levels = order_genes)

    ggplot2::ggplot(df, ggplot2::aes(x = .data$gene, y = .data$signed,
                                     fill = .data$alteration)) +
        ggplot2::geom_col() +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
        ggplot2::coord_flip() +
        ggplot2::scale_fill_manual(values = c(gain = "#B2182B",
                                              loss = "#2166AC")) +
        ggplot2::labs(x = NULL, y = "Fraction of samples", fill = NULL) +
        .oncoTheme()
}

#' Mutation frequency plot with confidence intervals
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param topN Number of genes to show.
#' @param nonsynonymousOnly Restrict to protein-altering variants.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' plotMutationFrequency(ose, topN = 8)
#'
#' @export
plotMutationFrequency <- function(object, topN = 20L,
                                  nonsynonymousOnly = TRUE) {
    df <- mutationFrequency(object, topN = topN,
                            nonsynonymousOnly = nonsynonymousOnly)
    df$gene <- factor(df$gene, levels = rev(df$gene))

    ggplot2::ggplot(df, ggplot2::aes(y = .data$gene, x = .data$frequency)) +
        ggplot2::geom_col(fill = "#4575B4") +
        ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ci_lower,
                                            xmax = .data$ci_upper),
                               orientation = "y", width = 0.25,
                               linewidth = 0.3) +
        ggplot2::labs(x = "Fraction of samples mutated", y = NULL) +
        .oncoTheme()
}
