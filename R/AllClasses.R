#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importClassesFrom S4Vectors DataFrame SimpleList
#' @importFrom S4Vectors DataFrame SimpleList metadata
#' @importFrom methods new validObject is as slot slot<- callNextMethod show
NULL

setClassUnion("SummarizedExperimentOrNULL",
    c("SummarizedExperiment", "NULL"))

#' Container for a cancer multi-omic study
#'
#' @description
#' `OncoSeqExperiment` extends \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} so that bulk
#' expression data, somatic mutations, gene-level copy-number alterations,
#' clinical follow-up and the results of every analysis performed on them live
#' in a single object. All layers are keyed on the same sample identifiers,
#' which are the column names of the object; validity checking enforces this.
#'
#' @slot mutations A \link[S4Vectors:DataFrame-class]{DataFrame} of somatic variants in long
#'   format, one row per variant per sample. Must contain at least the columns
#'   `sample` and `gene`; `variant_class` is used by the summarisation
#'   functions when present. Sample identifiers must be a subset of
#'   `colnames(object)`.
#' @slot cnv A \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} of gene-level copy-number
#'   calls (genes in rows, samples in columns), or `NULL`. The first assay is
#'   taken to hold discrete GISTIC-style calls in `-2:2`. Column names must be
#'   a subset of `colnames(object)`.
#' @slot analyses A \link[S4Vectors:SimpleList-class]{SimpleList} registry of analysis results.
#'   Each element is named and carries an `analysis` attribute describing how
#'   it was produced.
#' @slot provenance A list recording the import source and every processing
#'   step applied to the object.
#'
#' @section Accessors:
#' `mutations()`, `cnv()`, `analyses()` and `provenance()` return the
#' corresponding slot; each has a replacement counterpart. All
#' \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} accessors (`assay()`, `colData()`,
#' `rowData()`, `rowRanges()`) work unchanged.
#'
#' @section Subsetting:
#' `x[i, j]` subsets genes and samples. Subsetting samples also filters the
#' mutation table and the copy-number layer to the retained samples, so the
#' object stays internally consistent.
#'
#' @seealso [OncoSeqExperiment()] to construct an object, [importCounts()] and
#'   [importGDC()] to build one from data.
#'
#' @name OncoSeqExperiment-class
#' @aliases OncoSeqExperiment-class
#' @exportClass OncoSeqExperiment
setClass("OncoSeqExperiment",
    contains = "SummarizedExperiment",
    slots = c(
        mutations = "DataFrame",
        cnv = "SummarizedExperimentOrNULL",
        analyses = "SimpleList",
        provenance = "list"
    ),
    prototype = prototype(
        mutations = DataFrame(),
        cnv = NULL,
        analyses = SimpleList(),
        provenance = list()
    )
)

.validOncoSeqExperiment <- function(object) {
    msg <- character(0)
    samples <- colnames(object)

    mut <- object@mutations
    if (nrow(mut) > 0L) {
        required <- c("sample", "gene")
        missing_cols <- setdiff(required, colnames(mut))
        if (length(missing_cols) > 0L) {
            msg <- c(msg, paste0("'mutations' is missing required column(s): ",
                paste(missing_cols, collapse = ", ")))
        } else if (!is.null(samples)) {
            unknown <- setdiff(unique(as.character(mut$sample)), samples)
            if (length(unknown) > 0L) {
                msg <- c(msg, paste0(
                    "'mutations' refers to ", length(unknown),
                    " sample(s) absent from colnames(object), e.g. ",
                    paste(utils::head(unknown, 3L), collapse = ", ")))
            }
        }
    }

    cn <- object@cnv
    if (!is.null(cn)) {
        if (length(SummarizedExperiment::assays(cn)) == 0L) {
            msg <- c(msg, "'cnv' must contain at least one assay")
        }
        if (!is.null(samples)) {
            unknown <- setdiff(colnames(cn), samples)
            if (length(unknown) > 0L) {
                msg <- c(msg, paste0(
                    "'cnv' has ", length(unknown),
                    " column(s) absent from colnames(object), e.g. ",
                    paste(utils::head(unknown, 3L), collapse = ", ")))
            }
        }
    }

    an <- object@analyses
    if (length(an) > 0L) {
        nms <- names(an)
        if (is.null(nms) || any(!nzchar(nms))) {
            msg <- c(msg, "every element of 'analyses' must be named")
        } else if (anyDuplicated(nms) > 0L) {
            msg <- c(msg, "names of 'analyses' must be unique")
        }
    }

    if (length(msg) == 0L) TRUE else msg
}

setValidity("OncoSeqExperiment", .validOncoSeqExperiment)

#' Construct an OncoSeqExperiment
#'
#' @description
#' Builds an \linkS4class{OncoSeqExperiment} from an expression matrix or an
#' existing \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}, optionally attaching somatic
#' mutations, gene-level copy number and clinical data for the same samples.
#'
#' @param assays A named list of matrices with genes in rows and samples in
#'   columns, or a \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}. When a plain list is
#'   given the first element is treated as raw counts and should be named
#'   `counts`.
#' @param colData A `data.frame` or \link[S4Vectors:DataFrame-class]{DataFrame} of sample
#'   annotation whose rows correspond to the columns of `assays`. Ignored when
#'   `assays` is a \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}.
#' @param rowData Optional feature annotation.
#' @param mutations Optional somatic variant table in long format with columns
#'   `sample` and `gene` (see the class documentation).
#' @param cnv Optional gene-level copy-number data: a matrix or
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} with samples in columns.
#' @param provenance Optional named list describing the data source.
#' @param ... Further arguments passed to
#'   [SummarizedExperiment::SummarizedExperiment()].
#'
#' @return An \linkS4class{OncoSeqExperiment}.
#'
#' @examples
#' counts <- matrix(rpois(200, 50), nrow = 20,
#'                  dimnames = list(paste0("G", 1:20), paste0("S", 1:10)))
#' meta <- data.frame(group = rep(c("tumour", "normal"), each = 5),
#'                    row.names = colnames(counts))
#' muts <- data.frame(sample = c("S1", "S2"), gene = c("TP53", "RB1"),
#'                    variant_class = c("Missense_Mutation", "Nonsense_Mutation"))
#' ose <- OncoSeqExperiment(list(counts = counts), colData = meta,
#'                          mutations = muts)
#' ose
#'
#' @export
OncoSeqExperiment <- function(assays, colData = NULL, rowData = NULL,
                              mutations = NULL, cnv = NULL,
                              provenance = list(), ...) {
    if (is(assays, "SummarizedExperiment")) {
        se <- assays
    } else {
        if (is.matrix(assays) || is(assays, "Matrix")) {
            assays <- list(counts = assays)
        }
        if (!is.list(assays) || length(assays) == 0L) {
            stop("'assays' must be a matrix, a non-empty named list of ",
                 "matrices, or a SummarizedExperiment", call. = FALSE)
        }
        if (is.null(names(assays))) names(assays) <- "counts"
        args <- list(assays = assays)
        if (!is.null(colData)) args$colData <- .asDataFrame(colData)
        if (!is.null(rowData)) args$rowData <- .asDataFrame(rowData)
        se <- do.call(SummarizedExperiment::SummarizedExperiment,
                      c(args, list(...)))
    }

    if (is.null(colnames(se))) {
        colnames(se) <- paste0("sample", seq_len(ncol(se)))
    }
    if (is.null(rownames(se))) {
        stop("the expression data must have row names (feature identifiers)",
             call. = FALSE)
    }

    mut <- .normaliseMutationTable(mutations)
    cnvSE <- .asCNVSummarizedExperiment(cnv)

    out <- new("OncoSeqExperiment", se,
        mutations = mut,
        cnv = cnvSE,
        analyses = SimpleList(),
        provenance = c(list(created = Sys.time()), provenance))
    validObject(out)
    out
}
