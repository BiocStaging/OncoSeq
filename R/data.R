#' A small simulated cancer multi-omic dataset
#'
#' @description
#' Returns [simExample], a compact \linkS4class{OncoSeqExperiment} with
#' expression counts, somatic mutations, gene-level copy number and
#' overall-survival follow-up, optionally trimmed to fewer features or
#' samples. It is the object used by the examples and unit tests.
#'
#' @details
#' The data are simulated, not measured: a negative-binomial count model
#' with a planted differential signal between tumour and normal, mutations
#' planted at frequencies resembling those of urothelial carcinoma, and
#' survival times from an exponential model whose rate depends on `TP53`
#' mutation status. Gene symbols are real so that gene sets and driver
#' names behave sensibly in examples.
#'
#' The object is built once, by `data-raw/make-simExample.R`, and shipped
#' with the package. This function therefore draws no random numbers and
#' leaves the caller's random stream untouched; repeated calls return
#' identical objects.
#'
#' @param nGenes Number of features to keep, taken from the start of the
#'   object so that the named driver genes are always present. `NULL`
#'   keeps all 300.
#' @param nPerGroup Number of tumour and of normal samples to keep. `NULL`
#'   keeps all 20 of each.
#'
#' @return An \linkS4class{OncoSeqExperiment} with a `counts` assay,
#'   `condition`, `age`, `stage`, `os_time` and `os_event` in `colData()`,
#'   a mutation table and a copy-number layer.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose
#' table(SummarizedExperiment::colData(ose)$condition)
#' dim(exampleOncoSeq(nGenes = 50, nPerGroup = 5))
#'
#' @seealso [simExample] for the full object, [importCounts()] to build one
#'   from your own data.
#' @export
exampleOncoSeq <- function(nGenes = NULL, nPerGroup = NULL) {
    simExample <- NULL  # satisfy R CMD check; replaced by the loaded object
    utils::data("simExample", package = "OncoSeq", envir = environment())
    object <- simExample

    if (!is.null(nGenes)) {
        nGenes <- as.integer(nGenes)
        if (nGenes < 1L || nGenes > nrow(object)) {
            stop("'nGenes' must be between 1 and ", nrow(object),
                 call. = FALSE)
        }
        object <- object[seq_len(nGenes), ]
    }
    if (!is.null(nPerGroup)) {
        nPerGroup <- as.integer(nPerGroup)
        condition <- SummarizedExperiment::colData(object)$condition
        keep <- unlist(lapply(split(seq_len(ncol(object)), condition),
                              function(idx) utils::head(idx, nPerGroup)),
                       use.names = FALSE)
        if (length(keep) < 2L * nPerGroup) {
            stop("'nPerGroup' must be between 1 and ",
                 min(table(condition)), call. = FALSE)
        }
        object <- object[, sort(keep)]
    }
    object
}

#' Bladder urothelial carcinoma multi-omic example cohort
#'
#' @description
#' A real, trimmed TCGA bladder cancer cohort assembled from the cBioPortal
#' PanCancer Atlas study `blca_tcga_pan_can_atlas_2018` with
#' [importCBioPortal()]. It carries mRNA expression, somatic mutations,
#' discrete copy-number calls and clinical follow-up for the same 400
#' tumours, and is the dataset used throughout the vignette for the
#' genomic and survival sections.
#'
#' @details
#' Genes were restricted to a 420-gene panel: the 400 genes mutated in the
#' largest number of samples in the study, plus a curated list of
#' urothelial-carcinoma drivers and luminal/basal markers. Samples were
#' restricted to the 400 tumours carrying expression, mutation,
#' copy-number and overall-survival data simultaneously.
#'
#' The `expression` assay holds RSEM-normalised values, not raw counts, so
#' count-based engines will refuse to model it; use `engine = "limma"` with
#' `voom = FALSE`, or [exampleOncoSeq()] when a count-scale demonstration
#' is needed. Copy-number values are GISTIC calls in `-2:2`. Survival is
#' harmonised into `os_time` (months) and `os_event` (1 = death).
#'
#' @format An \linkS4class{OncoSeqExperiment} with 420 features and 400
#'   samples:
#' \describe{
#'   \item{assay `expression`}{RSEM-normalised mRNA abundance.}
#'   \item{`colData`}{clinical attributes including `AGE`, `SEX`,
#'     `AJCC_PATHOLOGIC_TUMOR_STAGE`, `SUBTYPE`, `TMB_NONSYNONYMOUS`,
#'     `ANEUPLOIDY_SCORE`, and the harmonised `os_time`/`os_event`.}
#'   \item{`mutations`}{14,337 somatic variants with `sample`, `gene`,
#'     `variant_class`, position and protein change.}
#'   \item{`cnv`}{discrete GISTIC gene-level copy-number calls.}
#' }
#'
#' @source cBioPortal study `blca_tcga_pan_can_atlas_2018`
#'   (<https://www.cbioportal.org/study/summary?id=blca_tcga_pan_can_atlas_2018>),
#'   derived from the TCGA PanCancer Atlas. The build script is
#'   `data-raw/make-blcaExample.R`.
#'
#' @examples
#' data(blcaExample)
#' blcaExample
#' head(mutationFrequency(blcaExample, topN = 10))
#'
#' @seealso [exampleOncoSeq()] for a simulated count-scale object.
"blcaExample"


#' Simulated cancer multi-omic cohort
#'
#' @description
#' The object returned by [exampleOncoSeq()]: 300 features by 40 samples,
#' 20 simulated tumours and 20 simulated normals, with counts, somatic
#' mutations, gene-level copy number and overall survival. The values are
#' simulated; see [exampleOncoSeq()] for the generative model and
#' `data-raw/make-simExample.R` for the code that produced it.
#'
#' @format An \linkS4class{OncoSeqExperiment} with 300 features and 40
#'   samples.
#'
#' @source Simulated by `data-raw/make-simExample.R`.
#'
#' @examples
#' data(simExample)
#' simExample
#'
#' @seealso [blcaExample] for a real cohort.
"simExample"
