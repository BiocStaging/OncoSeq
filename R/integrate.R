#' Mutation status of a gene as a factor
#'
#' @description
#' Returns a two-level factor giving, for every sample in the object,
#' whether the gene carries a somatic mutation. Samples absent from the
#' mutation table are treated as wild-type, which is correct when the whole
#' cohort was sequenced; check this assumption for cohorts with partial
#' sequencing.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param gene Gene symbol.
#' @param nonsynonymousOnly Count only protein-altering variants.
#'
#' @return A factor with levels `wild-type` and `mutant`, named by sample.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' table(mutationStatus(ose, "TP53"))
#'
#' @export
mutationStatus <- function(object, gene, nonsynonymousOnly = TRUE) {
    mut <- .mutationTable(object, nonsynonymousOnly)
    if (!gene %in% mut$gene) {
        warning("gene ", sQuote(gene), " carries no variant in this cohort; ",
                "all samples are wild-type", call. = FALSE)
    }
    hit <- unique(mut$sample[mut$gene == gene])
    stats::setNames(
        factor(ifelse(colnames(object) %in% hit, "mutant", "wild-type"),
               levels = c("wild-type", "mutant")),
        colnames(object))
}

#' Differential expression between mutant and wild-type tumours
#'
#' @description
#' Adds the mutation status of a driver gene to `colData()` and runs
#' [diffExpress()] comparing mutant with wild-type samples, optionally
#' adjusting for clinical covariates. This is the standard way to ask what
#' transcriptional programme accompanies an alteration.
#'
#' @param object An \linkS4class{OncoSeqExperiment} with expression and
#'   mutation layers.
#' @param gene Driver gene whose status defines the comparison.
#' @param covariates Optional `colData()` columns to adjust for.
#' @param engine Differential-expression engine; see [diffExpress()].
#' @param nonsynonymousOnly Count only protein-altering variants.
#' @param minGroupSize Refuse to run when either group is smaller than this.
#' @param name Name under which the result is stored.
#' @param ... Passed to [diffExpress()].
#'
#' @return The object with the result stored in the registry and a
#'   `<gene>_status` column added to `colData()`.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose <- filterLowExpression(ose, group = "condition", verbose = FALSE)
#' ose <- integrateMutationExpression(ose, "TP53", engine = "edgeR")
#' head(getAnalysis(ose, "de_TP53"))
#'
#' @seealso [diffExpress()], [mutationStatus()]
#' @export
integrateMutationExpression <- function(object, gene, covariates = NULL,
                                        engine = "edgeR",
                                        nonsynonymousOnly = TRUE,
                                        minGroupSize = 3L,
                                        name = paste0("de_", gene), ...) {
    status <- mutationStatus(object, gene, nonsynonymousOnly)
    sizes <- table(status)
    if (any(sizes < minGroupSize)) {
        stop("comparison on ", sQuote(gene), " has groups of size ",
             paste(paste0(names(sizes), "=", as.integer(sizes)),
                   collapse = ", "),
             "; at least ", minGroupSize, " per group are required",
             call. = FALSE)
    }
    col <- paste0(make.names(gene), "_status")
    SummarizedExperiment::colData(object)[[col]] <- status

    design <- stats::as.formula(
        paste("~", paste(c(covariates, col), collapse = " + ")))
    diffExpress(object, design,
                contrast = c(col, "mutant", "wild-type"),
                engine = engine, name = name, ...)
}

#' Restrict an object to samples present in every layer
#'
#' @description
#' Intersects the sample sets of the expression, mutation, copy-number and
#' clinical layers so that downstream multi-omic analyses operate on a
#' consistent cohort, and reports what was dropped.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param require Layers that a sample must have data for. Any of
#'   `"mutation"`, `"cnv"`, `"survival"`.
#' @param timeColumn,eventColumn Columns checked when `"survival"` is
#'   required.
#' @param verbose Report the number of samples retained.
#'
#' @return The subset \linkS4class{OncoSeqExperiment}.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' alignLayers(ose, require = c("mutation", "survival"))
#'
#' @export
alignLayers <- function(object, require = c("mutation", "cnv"),
                        timeColumn = "os_time", eventColumn = "os_event",
                        verbose = TRUE) {
    require <- match.arg(require, c("mutation", "cnv", "survival"),
                         several.ok = TRUE)
    keep <- colnames(object)
    cd <- SummarizedExperiment::colData(object)

    if ("mutation" %in% require) {
        mut <- mutations(object)
        if (nrow(mut) == 0L) {
            stop("'mutation' was required but the object has no mutation ",
                 "layer", call. = FALSE)
        }
        keep <- intersect(keep, unique(as.character(mut$sample)))
    }
    if ("cnv" %in% require) {
        cn <- cnv(object)
        if (is.null(cn)) {
            stop("'cnv' was required but the object has no copy-number layer",
                 call. = FALSE)
        }
        keep <- intersect(keep, colnames(cn))
    }
    if ("survival" %in% require) {
        for (nm in c(timeColumn, eventColumn)) {
            if (!nm %in% colnames(cd)) {
                stop("'survival' was required but colData has no column ",
                     sQuote(nm), call. = FALSE)
            }
        }
        complete <- !is.na(cd[[timeColumn]]) & !is.na(cd[[eventColumn]])
        keep <- intersect(keep, colnames(object)[complete])
    }

    if (length(keep) == 0L) {
        stop("no sample carries all of the required layers (",
             paste(require, collapse = ", "), ")", call. = FALSE)
    }
    if (verbose) {
        message("alignLayers: keeping ", length(keep), " of ", ncol(object),
                " samples with ", paste(require, collapse = " + "), " data")
    }
    out <- object[, keep]
    .addProvenance(out, "alignLayers", require = require,
                   kept = length(keep))
}
