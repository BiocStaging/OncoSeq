test_that("importCounts accepts matrices and matches sample annotation", {
    counts <- matrix(rpois(300, 40), nrow = 30,
        dimnames = list(paste0("G", 1:30), paste0("S", 1:10)))
    meta <- data.frame(condition = rep(c("a", "b"), each = 5),
                       row.names = colnames(counts))
    ose <- importCounts(counts, meta)
    expect_s4_class(ose, "OncoSeqExperiment")
    expect_identical(SummarizedExperiment::colData(ose)$condition,
                     meta$condition)
})

test_that("importCounts intersects mismatched sample sets", {
    counts <- matrix(rpois(300, 40), nrow = 30,
        dimnames = list(paste0("G", 1:30), paste0("S", 1:10)))
    meta <- data.frame(condition = rep(c("a", "b"), each = 4),
                       row.names = paste0("S", 3:10))
    ose <- importCounts(counts, meta)
    expect_identical(ncol(ose), 8L)
    expect_error(importCounts(counts, meta, intersectSamples = FALSE),
                 "intersectSamples")
})

test_that("importCounts rejects non-numeric input", {
    df <- data.frame(gene = paste0("G", 1:5), S1 = 1:5, S2 = 6:10)
    expect_error(importCounts(as.matrix(df)), "must be numeric")
})

test_that("importCounts reads delimited files", {
    counts <- matrix(rpois(60, 40), nrow = 6,
        dimnames = list(paste0("G", 1:6), paste0("S", 1:10)))
    f <- tempfile(fileext = ".tsv")
    df <- data.frame(gene = rownames(counts), counts, check.names = FALSE)
    utils::write.table(df, f, sep = "\t", quote = FALSE, row.names = FALSE)
    ose <- importCounts(f)
    expect_identical(dim(ose), c(6L, 10L))
    expect_identical(rownames(ose), rownames(counts))
})

test_that("importMAF harmonises a MAF file", {
    maf <- importMAF(system.file("extdata", "example_mutations.maf",
                                 package = "OncoSeq"))
    expect_s4_class(maf, "DFrame")
    expect_true(all(c("sample", "gene", "variant_class") %in% colnames(maf)))
    expect_gt(nrow(maf), 0L)
})

test_that("asOncoSeq coerces a SummarizedExperiment", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(rpois(100, 30), nrow = 10,
            dimnames = list(paste0("G", 1:10), paste0("S", 1:10)))))
    expect_s4_class(asOncoSeq(se), "OncoSeqExperiment")
    expect_error(asOncoSeq(1:10), "must be a SummarizedExperiment")
})

test_that("optional-dependency importers fail with an actionable message", {
    skip_if(requireNamespace("TCGAbiolinks", quietly = TRUE),
            "TCGAbiolinks is installed")
    expect_error(importGDC("TCGA-BLCA"), "BiocManager::install")
})
