test_that("survivalFit stratifies by mutation status", {
    fit <- survivalFit(ose_small, "TP53", by = "mutation")
    expect_s3_class(fit, "oncoseq_survival")
    expect_identical(fit$groups, c("wild-type", "mutant"))
    expect_true(fit$p_logrank >= 0 && fit$p_logrank <= 1)
    expect_output(print(fit), "Kaplan-Meier")
})

test_that("expression splits produce the expected group counts", {
    med <- survivalFit(ose_small, "MKI67", by = "expression",
                       split = "median")
    expect_identical(length(med$groups), 2L)
    expect_equal(sum(med$n), ncol(ose_small))
    ter <- survivalFit(ose_small, "MKI67", by = "expression",
                       split = "tertile")
    expect_identical(length(ter$groups), 3L)
    qua <- survivalFit(ose_small, "MKI67", by = "expression",
                       split = "quartile")
    expect_lt(sum(qua$n), ncol(ose_small))
})

test_that("the optimal cut-point is corrected for the search", {
    fit <- survivalFit(ose_small, "MKI67", by = "expression",
                       split = "optimal")
    expect_false(is.na(fit$cutpoint))
    expect_gte(fit$p_logrank_adjusted, fit$p_logrank)
})

test_that("missing survival columns give a clear error", {
    ose <- ose_small
    SummarizedExperiment::colData(ose)$os_time <- NULL
    expect_error(survivalFit(ose, "TP53", by = "mutation"),
                 "no column")
    expect_error(survivalFit(ose_small, "NOT_A_GENE"), "was not found")
})

test_that("coxModel accepts expression and mutation terms", {
    res <- coxModel(ose_small, c("expr:TP53", "mut:TP53", "age"),
                    name = NULL)
    expect_identical(nrow(res), 3L)
    expect_true(all(res$HR > 0))
    expect_true(all(res$ci_lower <= res$HR & res$HR <= res$ci_upper))
    expect_s3_class(attr(res, "model"), "coxph")
})

test_that("survivalScreen corrects across genes", {
    res <- survivalScreen(ose_small, genes = rownames(ose_small)[1:30],
                          name = NULL)
    expect_lte(nrow(res), 30L)
    expect_true(all(res$padj >= res$pvalue - 1e-12))
    expect_false(is.unsorted(res$padj))
})

test_that("survivalScreen stores its result when named", {
    out <- survivalScreen(ose_small, genes = rownames(ose_small)[1:20])
    expect_true("survival_screen" %in% getAnalysis(out))
})
