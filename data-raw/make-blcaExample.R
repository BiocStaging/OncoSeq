## Build data/blcaExample.rda from cBioPortal.
##
## Panel definition: the 400 genes mutated in the largest number of
## samples in the TCGA PanCancer Atlas bladder cohort, together with a
## curated list of urothelial-carcinoma drivers and luminal/basal markers.
## Expression (RSEM), discrete GISTIC copy number, somatic mutations and
## the clinical table are then fetched for that panel with the package's
## own importCBioPortal(), and the object is trimmed to the samples that
## carry all four layers.
##
## Run with:  Rscript data-raw/make-blcaExample.R

library(OncoSeq)
library(httr2)
library(jsonlite)
library(SummarizedExperiment)

STUDY <- "blca_tcga_pan_can_atlas_2018"
BASE <- "https://www.cbioportal.org/api"
N_PANEL <- 400L

post_json <- function(path, body) {
    req <- request(BASE)
    req <- req_url_path_append(req, path)
    req <- req_body_json(req, body, auto_unbox = FALSE)
    req <- req_retry(req, max_tries = 4L)
    req <- req_timeout(req, 900)
    fromJSON(resp_body_string(req_perform(req)), simplifyDataFrame = TRUE)
}

message("fetching the study mutation profile to define the gene panel")
muts <- post_json(
    paste0("molecular-profiles/", STUDY, "_mutations/mutations/fetch",
           "?projection=SUMMARY"),
    list(sampleListId = paste0(STUDY, "_sequenced")))

pairs <- unique(muts[, c("entrezGeneId", "sampleId")])
recurrence <- sort(table(pairs$entrezGeneId), decreasing = TRUE)
panelEntrez <- as.integer(names(recurrence)[seq_len(N_PANEL)])

symbols <- post_json("genes/fetch?geneIdType=ENTREZ_GENE_ID",
                     as.character(panelEntrez))
panel <- sort(unique(symbols$hugoGeneSymbol))

curated <- c("TP53", "RB1", "FGFR3", "ERBB2", "ERBB3", "PIK3CA", "KDM6A",
             "ARID1A", "KMT2D", "STAG2", "CDKN2A", "MDM2", "EGFR", "CCND1",
             "E2F3", "KRT5", "KRT6A", "KRT14", "KRT20", "PPARG", "GATA3",
             "FOXA1", "UPK1A", "UPK2", "MKI67", "CDH1", "VIM", "ESR2",
             "TSC1", "FGFR1", "MYC", "YAP1", "SHH", "EGFR")
panel <- sort(unique(c(panel, curated)))
message("panel size: ", length(panel), " genes")

message("importing the study through importCBioPortal()")
ose <- importCBioPortal(STUDY, genes = panel, chunkSize = 200L)

message("restricting to samples with expression, mutation, CNV and survival")
ose <- alignLayers(ose, require = c("mutation", "cnv", "survival"))

## Drop clinical columns that are entirely missing or constant, and keep a
## readable subset so that the shipped object stays small.
cd <- as.data.frame(colData(ose))
keepCols <- c("sampleId", "patientId", "AGE", "SEX", "RACE",
              "AJCC_PATHOLOGIC_TUMOR_STAGE", "GRADE",
              "SUBTYPE", "TMB_NONSYNONYMOUS", "ANEUPLOIDY_SCORE",
              "MSI_SCORE_MANTIS", "OS_STATUS", "OS_MONTHS",
              "DFS_STATUS", "DFS_MONTHS",
              "os_time", "os_event", "os_time_unit")
keepCols <- intersect(keepCols, colnames(cd))
colData(ose) <- S4Vectors::DataFrame(cd[, keepCols, drop = FALSE],
                                     row.names = rownames(cd))

## Round the expression values; RSEM estimates carry no meaningful
## precision beyond three decimals and rounding shrinks the stored object.
assay(ose, "expression") <- round(assay(ose, "expression"), 3)

provenance(ose) <- c(provenance(ose),
    list(panel = "400 most recurrently mutated genes plus curated drivers",
         built_by = "data-raw/make-blcaExample.R",
         built_on = Sys.Date()))

blcaExample <- ose
message("final object: ", nrow(blcaExample), " genes x ",
        ncol(blcaExample), " samples; ",
        nrow(mutations(blcaExample)), " variants")

save(blcaExample, file = "data/blcaExample.rda",
     compress = "xz", compression_level = 9)
message("written to data/blcaExample.rda (",
        round(file.size("data/blcaExample.rda") / 1024^2, 2), " MB)")
