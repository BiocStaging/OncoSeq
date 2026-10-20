test_that("mutationStatus labels every sample", {
    st <- mutationStatus(ose_small, "TP53")
    expect_length(st, ncol(ose_small))
    expect_identical(levels(st), c("wild-type", "mutant"))
    expect_warning(mutationStatus(ose_small, "NOT_A_GENE"), "no variant")
})

test_that("integrateMutationExpression runs a stratified comparison", {
    out <- integrateMutationExpression(ose_filt, "TP53", engine = "edgeR",
                                       verbose = FALSE)
    expect_true("de_TP53" %in% getAnalysis(out))
    expect_true("TP53_status" %in%
        colnames(SummarizedExperiment::colData(out)))
})

test_that("alignLayers intersects the requested layers", {
    out <- alignLayers(ose_small, require = c("mutation", "survival"),
                       verbose = FALSE)
    expect_lte(ncol(out), ncol(ose_small))
    expect_true(all(colnames(out) %in%
        unique(as.character(mutations(ose_small)$sample))))
})

test_that("plot functions return ggplot objects", {
    ose <- diffExpress(ose_filt, ~ condition,
        contrast = c("condition", "tumour", "normal"), engine = "edgeR",
        verbose = FALSE)
    ose <- runGSEA(ose, readGMT(gmt_path), minSize = 3L)
    fit <- survivalFit(ose, "TP53", by = "mutation")
    plots <- list(
        plotSamplePCA(ose, colourBy = "condition"),
        plotVolcano(ose),
        plotMAplot(ose),
        plotExpressionHeatmap(ose, nTopGenes = 15L),
        plotKM(fit),
        plotOncoprint(ose, topN = 5L),
        plotEnrichmentBar(ose, "gsea"),
        plotCNVFrequency(ose, topN = 5L),
        plotMutationFrequency(ose, topN = 5L))
    for (p in plots) expect_s3_class(p, "ggplot")
    built <- ggplot2::ggplot_build(plots[[2]])
    expect_gt(nrow(built$data[[1]]), 0L)
})

test_that("plotKM rejects anything but a survivalFit result", {
    expect_error(plotKM(list()), "survivalFit")
})
