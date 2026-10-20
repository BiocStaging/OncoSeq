# OncoSeq 0.99.0

First release prepared for Bioconductor submission. The package was
rewritten from the ground up; nothing of the 0.1.0 prototype other than
the name and the intent remains.

## New container

* `OncoSeqExperiment`, an S4 class extending `SummarizedExperiment`, holds
  expression, somatic mutations, gene-level copy number, clinical
  follow-up, an analysis registry and a provenance log for one cohort.
  Validity checking keeps sample identifiers concordant across layers and
  subsetting propagates to every layer.

## New data import

* `importGDC()` retrieves a TCGA project through `TCGAbiolinks`. It
  requests the current `"STAR - Counts"` workflow; the prototype's
  `"HTSeq - Counts"` request was retired by the GDC at Data Release 32 and
  no longer returns data.
* `importCBioPortal()`, `cbioStudies()` and `cbioProfiles()` query the
  cBioPortal REST API for curated studies.
* `importCounts()`, `importTximport()`, `importMAF()` and `asOncoSeq()`
  cover local files and existing objects.

## New analysis

* Quality control: `qcMetrics()`, `filterLowExpression()`,
  `detectOutliers()`.
* Normalisation: `normalizeCounts()` with CPM, TMM, upper-quartile, RLE,
  VST, rlog and log2-TPM.
* Differential expression: `diffExpress()` over DESeq2, edgeR and
  limma-voom behind one contrast interface, plus `compareEngines()`.
* Enrichment: `runGSEA()` (fgsea) and `runORA()` (hypergeometric), with
  `readGMT()` and `msigdbrSets()`.
* Genomics: `mutationSummary()`, `mutationFrequency()`,
  `tumourMutationBurden()`, `oncoplotData()`, `cnvSummary()`,
  `cnvFrequency()`, `cnvExpressionCorrelation()`.
* Survival: `survivalFit()` with median, tertile, quartile and
  search-corrected optimal cut-points, `coxModel()` and
  `survivalScreen()`.
* Integration: `mutationStatus()`, `integrateMutationExpression()`,
  `alignLayers()`.
* Visualisation: nine `plot*()` functions, all returning `ggplot` objects.
* `runOncoSeq()` chains the central steps in one call.

## New data

* `blcaExample`, a 420-gene by 400-sample TCGA bladder cancer cohort with
  all four layers, and `exampleOncoSeq()`, a simulated count-scale object
  for examples and tests.
