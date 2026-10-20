#' OncoSeq: end-to-end analysis of cancer transcriptomes and genomes
#'
#' @description
#' OncoSeq couples bulk RNA-sequencing counts with somatic mutations,
#' gene-level copy number and clinical follow-up for one set of tumours in a
#' single container, and provides the analyses that a cancer genomics study
#' needs from import to figure.
#'
#' @section Workflow:
#' \enumerate{
#'   \item Import with [importGDC()], [importCBioPortal()],
#'     [importCounts()] or [importTximport()].
#'   \item Inspect with [qcMetrics()] and [detectOutliers()], filter with
#'     [filterLowExpression()].
#'   \item Transform with [normalizeCounts()].
#'   \item Test with [diffExpress()] and, if you wish, [compareEngines()].
#'   \item Interpret with [runGSEA()] and [runORA()].
#'   \item Characterise the genome with [mutationFrequency()],
#'     [tumourMutationBurden()], [cnvSummary()] and
#'     [cnvExpressionCorrelation()].
#'   \item Relate to outcome with [survivalFit()], [coxModel()] and
#'     [survivalScreen()].
#'   \item Integrate with [integrateMutationExpression()].
#'   \item Draw with [plotSamplePCA()], [plotVolcano()], [plotKM()],
#'     [plotOncoprint()] and the other `plot*` functions.
#' }
#' [runOncoSeq()] chains the central steps in one call.
#'
#' @name OncoSeq-package
#' @aliases OncoSeq
#' @keywords internal
"_PACKAGE"
