# OncoSeq

End-to-end analysis of cancer transcriptomes and genomes in R.

OncoSeq keeps bulk RNA-seq counts, somatic mutations, gene-level copy
number and clinical follow-up for one set of tumours in a single
container, and provides the analyses a cancer genomics study needs from
import to figure. Every result is stored back in the object, together
with a provenance log, so a complete study travels as one file.

## Installation

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("babasaraki/OncoSeq")
```

The TCGA importer additionally needs `TCGAbiolinks`:

```r
BiocManager::install("TCGAbiolinks")
```

## Quick start

```r
library(OncoSeq)

# A curated cBioPortal study, or your own counts with importCounts()
ose <- importCBioPortal("blca_tcga_pan_can_atlas_2018",
                        genes = c("TP53", "RB1", "FGFR3", "ERBB2"))

# ... or the shipped example cohort
ose <- exampleOncoSeq()

ose <- qcMetrics(ose)
ose <- filterLowExpression(ose, group = "condition")
ose <- normalizeCounts(ose, method = "vst")
ose <- diffExpress(ose, ~ condition,
                   contrast = c("condition", "tumour", "normal"),
                   engine = "DESeq2")

plotVolcano(ose)
head(getAnalysis(ose, "de"))
```

Or run the whole chain at once:

```r
res <- runOncoSeq(ose, ~ condition,
                  contrast = c("condition", "tumour", "normal"),
                  geneSets = readGMT("hallmark.gmt"))
names(analyses(res))
```

## What is in the box

| Stage | Functions |
| --- | --- |
| Import | `importGDC()`, `importCBioPortal()`, `importCounts()`, `importTximport()`, `importMAF()`, `asOncoSeq()` |
| Quality control | `qcMetrics()`, `filterLowExpression()`, `detectOutliers()` |
| Normalisation | `normalizeCounts()` — CPM, TMM, upper-quartile, RLE, VST, rlog, log2-TPM |
| Differential expression | `diffExpress()` over DESeq2 / edgeR / limma-voom, `compareEngines()` |
| Enrichment | `runGSEA()`, `runORA()`, `readGMT()`, `msigdbrSets()` |
| Mutations | `mutationSummary()`, `mutationFrequency()`, `tumourMutationBurden()`, `oncoplotData()` |
| Copy number | `cnvSummary()`, `cnvFrequency()`, `cnvExpressionCorrelation()` |
| Survival | `survivalFit()`, `coxModel()`, `survivalScreen()` |
| Integration | `mutationStatus()`, `integrateMutationExpression()`, `alignLayers()` |
| Plots | `plotSamplePCA()`, `plotVolcano()`, `plotMAplot()`, `plotExpressionHeatmap()`, `plotKM()`, `plotOncoprint()`, `plotEnrichmentBar()`, `plotCNVFrequency()`, `plotMutationFrequency()` |
| Pipeline | `runOncoSeq()` |

## The container

`OncoSeqExperiment` extends `SummarizedExperiment`, so `assay()`,
`colData()`, `rowData()` and friends work as usual. It adds:

- `mutations(x)` — somatic variants in long format
- `cnv(x)` — gene-level copy-number calls
- `analyses(x)` — every result you have computed
- `provenance(x)` — what was done to the object, in order

Subsetting samples propagates to every layer, and validity checking keeps
the sample identifiers concordant.

## Data

- `exampleOncoSeq()` — a simulated count-scale cohort for examples and
  teaching
- `blcaExample` — a real 420-gene by 400-sample slice of the TCGA bladder
  urothelial carcinoma cohort, with expression, mutations, copy number and
  survival

## Documentation

```r
vignette("OncoSeq")
```

## Licence

GPL-3.
