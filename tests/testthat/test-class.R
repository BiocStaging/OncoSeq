test_that("the constructor builds a valid object", {
    counts <- matrix(rpois(200, 40), nrow = 20,
        dimnames = list(paste0("G", 1:20), paste0("S", 1:10)))
    ose <- OncoSeqExperiment(list(counts = counts))
    expect_s4_class(ose, "OncoSeqExperiment")
    expect_true(validObject(ose))
    expect_identical(dim(ose), c(20L, 10L))
    expect_identical(nrow(mutations(ose)), 0L)
    expect_null(cnv(ose))
})

test_that("row names are required and column names are generated", {
    counts <- matrix(rpois(100, 40), nrow = 10)
    expect_error(OncoSeqExperiment(list(counts = counts)), "row names")
    rownames(counts) <- paste0("G", 1:10)
    ose <- OncoSeqExperiment(list(counts = counts))
    expect_identical(colnames(ose), paste0("sample", 1:10))
})

test_that("validity rejects layers with unknown samples", {
    counts <- matrix(rpois(100, 40), nrow = 10,
        dimnames = list(paste0("G", 1:10), paste0("S", 1:10)))
    bad <- data.frame(sample = "NOT_A_SAMPLE", gene = "TP53")
    expect_error(OncoSeqExperiment(list(counts = counts), mutations = bad),
                 "absent from colnames")

    badCNV <- matrix(0, nrow = 2, ncol = 2,
        dimnames = list(c("TP53", "RB1"), c("X1", "X2")))
    ose <- OncoSeqExperiment(list(counts = counts))
    expect_error(cnv(ose) <- badCNV, "absent from colnames")
})

test_that("MAF-style column names are harmonised", {
    counts <- matrix(rpois(100, 40), nrow = 10,
        dimnames = list(paste0("G", 1:10), paste0("S", 1:10)))
    maf <- data.frame(Hugo_Symbol = c("TP53", "RB1"),
                      Tumor_Sample_Barcode = c("S1", "S2"),
                      Variant_Classification = "Missense_Mutation")
    ose <- OncoSeqExperiment(list(counts = counts), mutations = maf)
    expect_true(all(c("sample", "gene", "variant_class") %in%
                    colnames(mutations(ose))))
    expect_identical(as.character(mutations(ose)$gene), c("TP53", "RB1"))
})

test_that("a mutation table without recognisable columns is rejected", {
    counts <- matrix(rpois(100, 40), nrow = 10,
        dimnames = list(paste0("G", 1:10), paste0("S", 1:10)))
    expect_error(
        OncoSeqExperiment(list(counts = counts),
                          mutations = data.frame(a = 1, b = 2)),
        "needs column")
})

test_that("subsetting samples propagates to every layer", {
    keep <- colnames(ose_small)[1:8]
    sub <- ose_small[, keep]
    expect_identical(ncol(sub), 8L)
    expect_true(all(as.character(mutations(sub)$sample) %in% keep))
    expect_true(all(colnames(cnv(sub)) %in% keep))
    expect_true(validObject(sub))
})

test_that("the analysis registry stores and retrieves results", {
    expect_identical(getAnalysis(ose_small), character(0))
    expect_error(getAnalysis(ose_small, "nope"), "registry is empty")
    ose <- diffExpress(ose_filt, ~ condition,
        contrast = c("condition", "tumour", "normal"),
        engine = "edgeR", verbose = FALSE)
    expect_true("de" %in% getAnalysis(ose))
    expect_error(getAnalysis(ose, "nope"), "Available: de")
})

test_that("show prints the layers", {
    expect_output(show(ose_small), "class: OncoSeqExperiment")
    expect_output(show(ose_small), "mutations: ")
    expect_output(show(ose_small), "cnv: ")
})
