ose_small <- exampleOncoSeq(nGenes = 120L, nPerGroup = 12L)
ose_filt <- filterLowExpression(ose_small, group = "condition",
                                verbose = FALSE)
gmt_path <- system.file("extdata", "example_sets.gmt", package = "OncoSeq")
