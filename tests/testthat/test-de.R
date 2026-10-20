test_that("every engine returns the harmonised columns", {
    expected <- c("feature", "baseMean", "log2FC", "lfcSE", "stat",
                  "pvalue", "padj", "engine")
    for (e in c("DESeq2", "edgeR", "limma")) {
        ose <- diffExpress(ose_filt, ~ condition,
            contrast = c("condition", "tumour", "normal"),
            engine = e, name = "de", verbose = FALSE)
        res <- getAnalysis(ose, "de")
        expect_identical(colnames(res), expected, info = e)
        expect_identical(nrow(res), nrow(ose_filt), info = e)
        expect_identical(unique(res$engine), e)
        expect_false(is.unsorted(res$padj, na.rm = TRUE))
    }
})

test_that("the contrast direction is consistent across engines", {
    fcs <- lapply(c("DESeq2", "edgeR", "limma"), function(e) {
        res <- getAnalysis(diffExpress(ose_filt, ~ condition,
            contrast = c("condition", "tumour", "normal"),
            engine = e, verbose = FALSE), "de")
        stats::setNames(res$log2FC, res$feature)
    })
    common <- Reduce(intersect, lapply(fcs, names))
    expect_gt(stats::cor(fcs[[1]][common], fcs[[2]][common]), 0.9)
    expect_gt(stats::cor(fcs[[1]][common], fcs[[3]][common]), 0.9)
})

test_that("reversing the contrast flips the sign", {
    a <- getAnalysis(diffExpress(ose_filt, ~ condition,
        contrast = c("condition", "tumour", "normal"), engine = "edgeR",
        verbose = FALSE), "de")
    b <- getAnalysis(diffExpress(ose_filt, ~ condition,
        contrast = c("condition", "normal", "tumour"), engine = "edgeR",
        verbose = FALSE), "de")
    idx <- match(a$feature, b$feature)
    expect_equal(a$log2FC, -b$log2FC[idx], tolerance = 1e-8)
})

test_that("misspecified designs and contrasts are caught", {
    expect_error(diffExpress(ose_filt, ~ nope,
        contrast = c("nope", "a", "b")), "not found in colData")
    expect_error(diffExpress(ose_filt, ~ condition,
        contrast = c("condition", "tumour", "nonexistent")),
        "not found in")
    expect_error(diffExpress(ose_filt, ~ condition), "contrast")
})

test_that("DESeq2 refuses non-integer input", {
    ose <- normalizeCounts(ose_filt, "logcpm")
    expect_error(diffExpress(ose, ~ condition,
        contrast = c("condition", "tumour", "normal"),
        engine = "DESeq2", assay = "logcpm"), "integer counts")
})

test_that("limma runs on a transformed assay without voom", {
    ose <- normalizeCounts(ose_filt, "logcpm")
    out <- diffExpress(ose, ~ condition,
        contrast = c("condition", "tumour", "normal"),
        engine = "limma", assay = "logcpm", voom = FALSE, verbose = FALSE)
    expect_identical(nrow(getAnalysis(out, "de")), nrow(ose))
})

test_that("a covariate can be adjusted for", {
    out <- diffExpress(ose_filt, ~ stage + condition,
        contrast = c("condition", "tumour", "normal"),
        engine = "edgeR", verbose = FALSE)
    expect_identical(nrow(getAnalysis(out, "de")), nrow(ose_filt))
})

test_that("compareEngines reports an overlap", {
    cmp <- compareEngines(ose_filt, ~ condition,
        contrast = c("condition", "tumour", "normal"))
    expect_named(cmp, c("merged", "overlap", "consensus", "alpha"))
    expect_true("consensus" %in% names(cmp$overlap))
    expect_lte(cmp$overlap[["consensus"]], min(cmp$overlap[1:3]))
})
