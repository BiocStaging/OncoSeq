test_that("readGMT parses the shipped collection", {
    sets <- readGMT(gmt_path)
    expect_type(sets, "list")
    expect_true(all(lengths(sets) > 0L))
    expect_true("CELL_CYCLE_CONTROL" %in% names(sets))
    expect_length(readGMT(gmt_path, minSize = 1000L), 0L)
})

test_that("runGSEA stores a harmonised table", {
    ose <- diffExpress(ose_filt, ~ condition,
        contrast = c("condition", "tumour", "normal"), engine = "edgeR",
        verbose = FALSE)
    out <- runGSEA(ose, readGMT(gmt_path), minSize = 3L)
    res <- getAnalysis(out, "gsea")
    expect_true(all(c("set", "size", "ES", "NES", "pvalue", "padj",
                      "leadingEdge", "method") %in% colnames(res)))
    expect_true(all(res$padj >= 0 & res$padj <= 1, na.rm = TRUE))
})

test_that("runORA computes the hypergeometric test correctly", {
    ose <- diffExpress(ose_filt, ~ condition,
        contrast = c("condition", "tumour", "normal"), engine = "edgeR",
        verbose = FALSE)
    out <- runORA(ose, readGMT(gmt_path), alpha = 0.5, minSize = 3L)
    res <- getAnalysis(out, "ora")
    de <- getAnalysis(ose, "de")
    k <- sum(de$padj < 0.5, na.rm = TRUE)
    N <- nrow(de)
    row <- res[1L, ]
    expect_equal(row$pvalue,
        stats::phyper(row$overlap - 1L, row$size, N - row$size, k,
                      lower.tail = FALSE))
    expect_equal(row$expected, k * row$size / N)
})

test_that("mismatched identifiers give an informative error", {
    ose <- diffExpress(ose_filt, ~ condition,
        contrast = c("condition", "tumour", "normal"), engine = "edgeR",
        verbose = FALSE)
    expect_error(runGSEA(ose, list(setA = c("ENSG000001", "ENSG000002"))),
                 "common namespace")
})
