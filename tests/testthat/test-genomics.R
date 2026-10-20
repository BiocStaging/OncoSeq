test_that("mutationSummary counts variants, genes and samples", {
    s <- mutationSummary(ose_small)
    expect_named(s, c("perSample", "perGene", "byClass"))
    expect_identical(nrow(s$perSample), ncol(ose_small))
    mut <- as.data.frame(mutations(ose_small))
    mut <- mut[mut$variant_class %in% OncoSeq:::NONSYNONYMOUS_CLASSES, ]
    expect_identical(sum(s$perSample$n_variants), nrow(mut))
})

test_that("mutationFrequency uses the full cohort as denominator", {
    f <- mutationFrequency(ose_small, topN = Inf)
    expect_true(all(f$n_samples == ncol(ose_small)))
    expect_equal(f$frequency, f$n_mutated / ncol(ose_small))
    expect_true(all(f$ci_lower <= f$frequency & f$frequency <= f$ci_upper))
    expect_false(is.unsorted(rev(f$frequency)))
})

test_that("tumourMutationBurden scales by capture size", {
    a <- tumourMutationBurden(ose_small, captureSize = 38, store = FALSE)
    b <- tumourMutationBurden(ose_small, captureSize = 19, store = FALSE)
    expect_equal(b$tmb, a$tmb * 2)
    expect_error(tumourMutationBurden(ose_small, captureSize = 0),
                 "positive number")
})

test_that("oncoplotData has the right shape and marks multi-hits", {
    m <- oncoplotData(ose_small, topN = 5L)
    expect_identical(ncol(m), ncol(ose_small))
    expect_lte(nrow(m), 5L)
    expect_true(all(m[m != ""] %in%
        c(OncoSeq:::NONSYNONYMOUS_CLASSES, "Multi_Hit")))
})

test_that("mutation functions error clearly without a mutation layer", {
    bare <- OncoSeqExperiment(list(counts =
        matrix(rpois(100, 30), nrow = 10,
               dimnames = list(paste0("G", 1:10), paste0("S", 1:10)))))
    expect_error(mutationSummary(bare), "no mutation data")
})

test_that("cnvSummary counts gains and losses", {
    cs <- cnvSummary(ose_small)
    mat <- SummarizedExperiment::assay(cnv(ose_small), 1L)
    expect_equal(sum(cs$perGene$n_gain), sum(mat >= 2))
    expect_equal(sum(cs$perGene$n_loss), sum(mat <= -2))
    expect_true(all(cs$perSample$fraction_altered <= 1))
})

test_that("cnvFrequency returns long format", {
    f <- cnvFrequency(ose_small, topN = 5L)
    expect_identical(sort(unique(f$alteration)), c("gain", "loss"))
    expect_true(all(f$frequency >= 0 & f$frequency <= 1))
})

test_that("cnvExpressionCorrelation tests shared genes only", {
    res <- cnvExpressionCorrelation(ose_small)
    shared <- intersect(rownames(ose_small), rownames(cnv(ose_small)))
    expect_true(all(res$gene %in% shared))
    expect_true(all(res$padj >= res$pvalue - 1e-12, na.rm = TRUE))
})

test_that("copy-number functions error clearly without a CNV layer", {
    bare <- OncoSeqExperiment(list(counts =
        matrix(rpois(100, 30), nrow = 10,
               dimnames = list(paste0("G", 1:10), paste0("S", 1:10)))))
    expect_error(cnvSummary(bare), "no copy-number data")
})
