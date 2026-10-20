test_that("qcMetrics adds the expected columns", {
    ose <- qcMetrics(ose_small)
    cd <- SummarizedExperiment::colData(ose)
    expect_true(all(c("library_size", "detected_genes",
                      "median_expression") %in% colnames(cd)))
    expect_equal(cd$library_size,
                 colSums(SummarizedExperiment::assay(ose_small, "counts")),
                 ignore_attr = TRUE)
    metrics <- qcMetrics(ose_small, store = FALSE)
    expect_s4_class(metrics, "DFrame")
})

test_that("filterLowExpression removes features and records provenance", {
    out <- filterLowExpression(ose_small, group = "condition",
                               verbose = FALSE)
    expect_lte(nrow(out), nrow(ose_small))
    expect_true(any(vapply(provenance(out)$steps, function(s) {
        identical(s$step, "filterLowExpression")
    }, logical(1))))
})

test_that("filterLowExpression errors on an unknown grouping column", {
    expect_error(filterLowExpression(ose_small, group = "nope"),
                 "not a column of colData")
})

test_that("detectOutliers returns one row per sample", {
    out <- detectOutliers(ose_small)
    expect_identical(nrow(out), ncol(ose_small))
    expect_true(is.logical(out$outlier))
})

test_that("each normalisation method adds a named assay", {
    for (m in c("cpm", "logcpm", "tmm", "upperquartile", "rle", "vst")) {
        out <- normalizeCounts(ose_filt, method = m)
        expect_true(m %in% SummarizedExperiment::assayNames(out),
                    info = m)
        expect_identical(dim(SummarizedExperiment::assay(out, m)),
                         dim(ose_filt))
    }
})

test_that("count-based normalisation rejects negative values", {
    ose <- ose_filt
    SummarizedExperiment::assay(ose, "counts")[1, 1] <- -5
    expect_error(normalizeCounts(ose, "cpm"), "non-negative counts")
})

test_that("log2tpm needs feature lengths", {
    expect_error(normalizeCounts(ose_filt, "log2tpm"), "feature lengths")
    ose <- ose_filt
    SummarizedExperiment::rowData(ose)$length <-
        sample(500:5000, nrow(ose), replace = TRUE)
    out <- normalizeCounts(ose, "log2tpm")
    expect_true("log2tpm" %in% SummarizedExperiment::assayNames(out))
})
