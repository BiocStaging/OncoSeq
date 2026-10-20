#' Run the full OncoSeq workflow in one call
#'
#' @description
#' Executes the standard sequence of steps on an existing object: quality
#' metrics, low-expression filtering, normalisation, differential
#' expression, enrichment, mutation burden and a survival screen, storing
#' every result in the analysis registry. Steps whose prerequisites are
#' absent are skipped with a message rather than failing, so the same call
#' works on expression-only and full multi-omic objects.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param design A one-sided formula for the differential-expression model.
#' @param contrast Contrast specification; see [diffExpress()].
#' @param geneSets Optional named list of gene sets enabling the enrichment
#'   step.
#' @param engine Differential-expression engine.
#' @param normalisation Normalisation method passed to [normalizeCounts()].
#' @param survivalGenes Genes to carry through [survivalScreen()]; `NULL`
#'   screens every feature, which is slow on a full transcriptome.
#' @param timeColumn,eventColumn Survival columns in `colData()`.
#' @param alpha Significance threshold used throughout.
#' @param verbose Report progress.
#'
#' @return The \linkS4class{OncoSeqExperiment} with results stored under the
#'   names `de`, `gsea`, `ora` and `survival_screen`, QC columns in
#'   `colData()` and the normalised assay added.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' sets <- readGMT(system.file("extdata", "example_sets.gmt",
#'                             package = "OncoSeq"))
#' res <- runOncoSeq(ose, ~ condition,
#'                   contrast = c("condition", "tumour", "normal"),
#'                   geneSets = sets, engine = "edgeR", verbose = FALSE)
#' names(analyses(res))
#'
#' @seealso [diffExpress()], [runGSEA()], [survivalScreen()]
#' @export
runOncoSeq <- function(object, design, contrast = NULL, geneSets = NULL,
                       engine = c("DESeq2", "edgeR", "limma"),
                       normalisation = "vst", survivalGenes = NULL,
                       timeColumn = "os_time", eventColumn = "os_event",
                       alpha = 0.05, verbose = TRUE) {
    stopifnot(is(object, "OncoSeqExperiment"))
    engine <- match.arg(engine)
    say <- function(...) if (verbose) message("[OncoSeq] ", ...)

    say("quality control")
    object <- qcMetrics(object)

    groupVar <- all.vars(stats::as.formula(design))
    groupVar <- if (length(groupVar) > 0L) utils::tail(groupVar, 1L) else NULL
    object <- filterLowExpression(object, group = groupVar,
                                  verbose = verbose)

    say("normalisation (", normalisation, ")")
    object <- try(normalizeCounts(object, method = normalisation),
                  silent = TRUE)
    if (inherits(object, "try-error")) {
        stop("normalisation failed: ", attr(object, "condition")$message,
             call. = FALSE)
    }

    say("differential expression (", engine, ")")
    object <- diffExpress(object, design, contrast = contrast,
                          engine = engine, alpha = alpha, verbose = verbose)

    if (!is.null(geneSets)) {
        say("gene set enrichment")
        object <- runGSEA(object, geneSets)
        ora <- try(runORA(object, geneSets, alpha = alpha), silent = TRUE)
        if (inherits(ora, "try-error")) {
            say("over-representation skipped: ",
                conditionMessage(attr(ora, "condition")))
        } else {
            object <- ora
        }
    } else {
        say("no gene sets supplied; enrichment skipped")
    }

    if (nrow(mutations(object)) > 0L) {
        say("tumour mutation burden")
        object <- tumourMutationBurden(object)
    }

    cd <- SummarizedExperiment::colData(object)
    if (all(c(timeColumn, eventColumn) %in% colnames(cd))) {
        say("survival screen")
        screen <- try(survivalScreen(object, genes = survivalGenes,
            timeColumn = timeColumn, eventColumn = eventColumn),
            silent = TRUE)
        if (inherits(screen, "try-error")) {
            say("survival screen skipped: ",
                conditionMessage(attr(screen, "condition")))
        } else {
            object <- screen
        }
    } else {
        say("no survival columns; survival screen skipped")
    }

    say("done; results: ", paste(names(analyses(object)), collapse = ", "))
    .addProvenance(object, "runOncoSeq", engine = engine,
                   normalisation = normalisation)
}
