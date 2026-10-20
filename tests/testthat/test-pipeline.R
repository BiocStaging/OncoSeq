test_that("runOncoSeq populates the registry", {
    out <- runOncoSeq(ose_small, ~ condition,
        contrast = c("condition", "tumour", "normal"),
        geneSets = readGMT(gmt_path), engine = "edgeR",
        normalisation = "logcpm",
        survivalGenes = rownames(ose_small)[1:20], verbose = FALSE)
    expect_s4_class(out, "OncoSeqExperiment")
    expect_true(all(c("de", "gsea") %in% getAnalysis(out)))
    expect_true("logcpm" %in% SummarizedExperiment::assayNames(out))
    expect_true("library_size" %in%
        colnames(SummarizedExperiment::colData(out)))
})

test_that("runOncoSeq skips steps whose inputs are absent", {
    bare <- OncoSeqExperiment(
        list(counts = matrix(rpois(600, 60), nrow = 30,
            dimnames = list(paste0("G", 1:30), paste0("S", 1:20)))),
        colData = data.frame(condition = rep(c("a", "b"), each = 10),
                             row.names = paste0("S", 1:20)))
    out <- runOncoSeq(bare, ~ condition, contrast = c("condition", "a", "b"),
        engine = "edgeR", normalisation = "logcpm", verbose = FALSE)
    expect_true("de" %in% getAnalysis(out))
    expect_false("gsea" %in% getAnalysis(out))
    expect_false("survival_screen" %in% getAnalysis(out))
})

test_that("provenance accumulates every step", {
    out <- runOncoSeq(ose_small, ~ condition,
        contrast = c("condition", "tumour", "normal"), engine = "edgeR",
        normalisation = "logcpm", verbose = FALSE)
    steps <- vapply(provenance(out)$steps, `[[`, character(1), "step")
    expect_true(all(c("qcMetrics", "filterLowExpression",
                      "normalizeCounts", "diffExpress") %in% steps))
})
