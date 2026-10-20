#' @importFrom stats model.matrix formula terms as.formula relevel
NULL

#' Differential expression with DESeq2, edgeR or limma
#'
#' @description
#' Fits a differential-expression model with one of three engines behind a
#' single interface and returns one harmonised result table. The engines are
#' the standard implementations: [DESeq2::DESeq()] with a negative-binomial
#' Wald test, [edgeR::glmQLFTest()] quasi-likelihood F tests, and
#' [limma::eBayes()] on [limma::voom()] precision weights (or on an already
#' transformed assay with `limma-trend`).
#'
#' @details
#' The contrast may be given in three ways:
#' \itemize{
#'   \item a character triple `c(variable, numerator, denominator)`, which
#'     tests `numerator - denominator` for that `colData()` variable;
#'   \item a coefficient name present in the design matrix, through `coef`;
#'   \item a numeric vector over the columns of the design matrix.
#' }
#' Additional terms in `design` are adjusted for. The same contrast
#' specification gives the same comparison in all three engines, so results
#' can be compared directly.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param design A one-sided formula over `colData()` columns, for example
#'   `~ condition` or `~ batch + condition`.
#' @param contrast A character triple `c(variable, numerator, denominator)`
#'   or a numeric contrast vector over the design-matrix columns.
#' @param coef Name or index of a single design-matrix coefficient to test;
#'   an alternative to `contrast`.
#' @param engine One of `"DESeq2"`, `"edgeR"` or `"limma"`.
#' @param assay Assay to model. Count engines default to `counts`; `limma`
#'   uses `counts` with voom unless `voom = FALSE`.
#' @param voom For `engine = "limma"`, apply [limma::voom()] to counts
#'   (`TRUE`), or fit `limma-trend` to the supplied transformed assay
#'   (`FALSE`). `NULL` chooses automatically from the assay.
#' @param shrink For `engine = "DESeq2"`, apply [DESeq2::lfcShrink()] with
#'   `apeglm`; requires the `apeglm` package and a `coef`.
#' @param name Name under which the result is stored in the analysis
#'   registry.
#' @param alpha Significance level used when reporting the number of
#'   differentially expressed features.
#' @param verbose Print a one-line summary of the fit.
#' @param ... Further arguments passed to the engine's fitting function
#'   ([DESeq2::DESeq()], [edgeR::glmQLFit()] or [limma::eBayes()]).
#'
#' @return The \linkS4class{OncoSeqExperiment} with the result table stored
#'   in `analyses(object)[[name]]`. Retrieve it with [getAnalysis()]. The
#'   table is a `data.frame` with columns `feature`, `baseMean`, `log2FC`,
#'   `lfcSE`, `stat`, `pvalue`, `padj` and `engine`, ordered by `padj`.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose <- filterLowExpression(ose, group = "condition", verbose = FALSE)
#' ose <- diffExpress(ose, ~ condition,
#'                    contrast = c("condition", "tumour", "normal"),
#'                    engine = "edgeR", name = "de_edger")
#' head(getAnalysis(ose, "de_edger"))
#'
#' @seealso [compareEngines()], [runGSEA()], [plotVolcano()]
#' @export
diffExpress <- function(object, design, contrast = NULL, coef = NULL,
                        engine = c("DESeq2", "edgeR", "limma"),
                        assay = NULL, voom = NULL, shrink = FALSE,
                        name = "de", alpha = 0.05, verbose = TRUE, ...) {
    stopifnot(is(object, "OncoSeqExperiment"))
    engine <- match.arg(engine)
    if (is.null(contrast) && is.null(coef)) {
        stop("supply either 'contrast' or 'coef' to define the comparison",
             call. = FALSE)
    }
    design <- stats::as.formula(design)
    cd <- as.data.frame(SummarizedExperiment::colData(object))
    .checkDesignVariables(design, cd)
    cd <- .dropUnusedLevels(cd, design)
    .checkContrastLevels(contrast, cd)

    mat <- as.matrix(.getAssay(object, assay))
    if (anyNA(mat)) {
        stop("the assay contains missing values, which none of the ",
             "differential-expression engines accept. Filter or impute ",
             "them first.", call. = FALSE)
    }

    res <- switch(engine,
        DESeq2 = .deDESeq2(mat, cd, design, contrast, coef, shrink, ...),
        edgeR = .deEdgeR(mat, cd, design, contrast, coef, ...),
        limma = .deLimma(mat, cd, design, contrast, coef, voom,
                         isCounts = .assayIsCounts(object, assay), ...))

    res$engine <- engine
    res <- res[order(res$padj, res$pvalue, na.last = TRUE), , drop = FALSE]
    rownames(res) <- NULL

    nSig <- sum(res$padj < alpha, na.rm = TRUE)
    if (verbose) {
        message(engine, ": ", nSig, " of ", nrow(res),
                " features with padj < ", alpha)
    }

    object <- .addAnalysis(object, name, res, analysis = "differential expression",
        engine = engine, design = deparse1(design), contrast = contrast,
        coef = coef, alpha = alpha, n_significant = nSig)
    .addProvenance(object, "diffExpress", engine = engine, name = name)
}

## Validate a c(variable, numerator, denominator) contrast before the
## engine sees it, so that all three engines give the same error.
.checkContrastLevels <- function(contrast, cd) {
    if (is.null(contrast) || is.numeric(contrast)) return(invisible(TRUE))
    if (length(contrast) != 3L) {
        stop("'contrast' must be c(variable, numerator, denominator) or a ",
             "numeric vector", call. = FALSE)
    }
    v <- as.character(contrast)[1L]
    if (!v %in% colnames(cd)) {
        stop("contrast variable ", sQuote(v), " is not a colData column; ",
             "available: ", paste(colnames(cd), collapse = ", "),
             call. = FALSE)
    }
    levs <- levels(factor(cd[[v]]))
    unknown <- setdiff(as.character(contrast)[2:3], levs)
    if (length(unknown) > 0L) {
        stop("level(s) ", paste(sQuote(unknown), collapse = ", "),
             " not found in ", sQuote(v), "; available: ",
             paste(levs, collapse = ", "), call. = FALSE)
    }
    invisible(TRUE)
}

.checkDesignVariables <- function(design, cd) {
    vars <- all.vars(design)
    missing_vars <- setdiff(vars, colnames(cd))
    if (length(missing_vars) > 0L) {
        stop("design variable(s) not found in colData(): ",
             paste(missing_vars, collapse = ", "), ". Available: ",
             paste(colnames(cd), collapse = ", "), call. = FALSE)
    }
    for (v in vars) {
        if (anyNA(cd[[v]])) {
            stop("design variable ", sQuote(v), " contains missing values ",
                 "in ", sum(is.na(cd[[v]])), " sample(s); subset the object ",
                 "to complete cases first", call. = FALSE)
        }
    }
    invisible(TRUE)
}

.dropUnusedLevels <- function(cd, design) {
    for (v in all.vars(design)) {
        if (is.character(cd[[v]])) cd[[v]] <- factor(cd[[v]])
        if (is.factor(cd[[v]])) cd[[v]] <- droplevels(cd[[v]])
    }
    cd
}

.assayIsCounts <- function(object, assay) {
    nms <- SummarizedExperiment::assayNames(object)
    used <- if (is.null(assay)) {
        if ("counts" %in% nms) "counts" else nms[1L]
    } else assay
    used == "counts"
}

## Build a numeric contrast over the columns of a design matrix from a
## c(variable, numerator, denominator) specification. With treatment
## contrasts the mean of a level is intercept + beta_level, so the
## difference of two levels is beta_num - beta_den.
.contrastVector <- function(dm, contrast, cd) {
    if (is.numeric(contrast)) {
        if (length(contrast) != ncol(dm)) {
            stop("a numeric 'contrast' must have one entry per design ",
                 "matrix column (", ncol(dm), ")", call. = FALSE)
        }
        return(contrast)
    }
    if (length(contrast) != 3L) {
        stop("'contrast' must be c(variable, numerator, denominator) or a ",
             "numeric vector", call. = FALSE)
    }
    v <- contrast[1L]
    if (!v %in% colnames(cd)) {
        stop("contrast variable ", sQuote(v), " is not a colData column",
             call. = FALSE)
    }
    levs <- levels(factor(cd[[v]]))
    for (lv in contrast[2:3]) {
        if (!lv %in% levs) {
            stop("level ", sQuote(lv), " not found in ", sQuote(v),
                 "; available: ", paste(levs, collapse = ", "),
                 call. = FALSE)
        }
    }
    cv <- rep(0, ncol(dm))
    names(cv) <- colnames(dm)
    numCol <- paste0(v, contrast[2L])
    denCol <- paste0(v, contrast[3L])
    if (numCol %in% colnames(dm)) cv[numCol] <- 1
    if (denCol %in% colnames(dm)) cv[denCol] <- -1
    if (all(cv == 0)) {
        stop("neither ", sQuote(numCol), " nor ", sQuote(denCol),
             " appears in the design matrix; check that 'design' includes ",
             sQuote(v), call. = FALSE)
    }
    cv
}

.resolveCoef <- function(dm, coef) {
    if (is.numeric(coef)) {
        if (coef < 1 || coef > ncol(dm)) {
            stop("'coef' index out of range", call. = FALSE)
        }
        return(as.integer(coef))
    }
    hit <- match(coef, colnames(dm))
    if (is.na(hit)) {
        stop("coefficient ", sQuote(coef), " is not a design matrix column; ",
             "available: ", paste(colnames(dm), collapse = ", "),
             call. = FALSE)
    }
    hit
}

.harmonise <- function(feature, baseMean, log2FC, lfcSE, stat, pvalue,
                       padj) {
    data.frame(
        feature = as.character(feature),
        baseMean = as.numeric(baseMean),
        log2FC = as.numeric(log2FC),
        lfcSE = as.numeric(lfcSE),
        stat = as.numeric(stat),
        pvalue = as.numeric(pvalue),
        padj = as.numeric(padj),
        stringsAsFactors = FALSE)
}

.deDESeq2 <- function(mat, cd, design, contrast, coef, shrink, ...) {
    counts <- mat
    if (!isTRUE(all.equal(counts, round(counts), tolerance = 1e-8))) {
        stop("DESeq2 requires integer counts; use engine = 'limma' for ",
             "continuous data", call. = FALSE)
    }
    storage.mode(counts) <- "integer"
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData = cd,
                                          design = design)
    dds <- DESeq2::DESeq(dds, quiet = TRUE, ...)

    if (!is.null(contrast) && !is.numeric(contrast)) {
        res <- DESeq2::results(dds, contrast = as.character(contrast))
    } else if (!is.null(contrast)) {
        res <- DESeq2::results(dds, contrast = contrast)
    } else {
        nm <- if (is.numeric(coef)) DESeq2::resultsNames(dds)[coef] else coef
        if (!nm %in% DESeq2::resultsNames(dds)) {
            stop("coefficient ", sQuote(nm), " is not among resultsNames(): ",
                 paste(DESeq2::resultsNames(dds), collapse = ", "),
                 call. = FALSE)
        }
        res <- DESeq2::results(dds, name = nm)
        if (isTRUE(shrink)) {
            .requireNamespaces("apeglm", "lfcShrink with type = 'apeglm'")
            res <- DESeq2::lfcShrink(dds, coef = nm, type = "apeglm",
                                     res = res, quiet = TRUE)
        }
    }
    res <- as.data.frame(res)
    .harmonise(rownames(res), res$baseMean, res$log2FoldChange, res$lfcSE,
               res$stat %||% rep(NA_real_, nrow(res)), res$pvalue, res$padj)
}

.deEdgeR <- function(mat, cd, design, contrast, coef, ...) {
    dm <- stats::model.matrix(design, data = cd)
    dge <- edgeR::DGEList(counts = mat)
    dge <- edgeR::normLibSizes(dge)
    dge <- edgeR::estimateDisp(dge, dm, robust = TRUE)
    fit <- edgeR::glmQLFit(dge, dm, robust = TRUE, ...)

    if (!is.null(contrast)) {
        cv <- .contrastVector(dm, contrast, cd)
        test <- edgeR::glmQLFTest(fit, contrast = cv)
    } else {
        test <- edgeR::glmQLFTest(fit, coef = .resolveCoef(dm, coef))
    }
    tt <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
    .harmonise(rownames(tt), 2^tt$logCPM, tt$logFC,
               rep(NA_real_, nrow(tt)), tt$F, tt$PValue, tt$FDR)
}

.deLimma <- function(mat, cd, design, contrast, coef, voom, isCounts, ...) {
    dm <- stats::model.matrix(design, data = cd)
    useVoom <- if (is.null(voom)) isCounts else isTRUE(voom)

    if (useVoom) {
        dge <- edgeR::DGEList(counts = mat)
        dge <- edgeR::normLibSizes(dge)
        v <- limma::voom(dge, dm)
        fit <- limma::lmFit(v, dm)
        trend <- FALSE
        baseMean <- 2^rowMeans(v$E)
    } else {
        fit <- limma::lmFit(mat, dm)
        trend <- TRUE
        baseMean <- rowMeans(mat)
    }

    if (!is.null(contrast)) {
        cv <- .contrastVector(dm, contrast, cd)
        fit <- limma::contrasts.fit(fit, contrasts = cv)
        fit <- limma::eBayes(fit, trend = trend, ...)
        tt <- limma::topTable(fit, number = Inf, sort.by = "none")
    } else {
        fit <- limma::eBayes(fit, trend = trend, ...)
        tt <- limma::topTable(fit, coef = .resolveCoef(dm, coef),
                              number = Inf, sort.by = "none")
    }
    .harmonise(rownames(tt), baseMean[rownames(tt)], tt$logFC,
               rep(NA_real_, nrow(tt)), tt$t, tt$P.Value, tt$adj.P.Val)
}

#' Run the same contrast through every engine and compare
#'
#' @description
#' Fits the same comparison with DESeq2, edgeR and limma-voom and returns the
#' overlap of their significant calls together with a merged table of
#' log fold changes. Concordance across engines is a practical check that a
#' signature is not an artefact of one model's assumptions.
#'
#' @inheritParams diffExpress
#' @param engines Engines to run.
#' @param alpha Adjusted-p-value threshold used for the overlap counts.
#'
#' @return A list with `merged` (one row per feature, log2FC and padj from
#'   each engine), `overlap` (a named integer vector of significant counts
#'   per engine and the size of their intersection) and `consensus` (the
#'   features significant in every engine).
#'
#' @examples
#' ose <- exampleOncoSeq()
#' ose <- filterLowExpression(ose, group = "condition", verbose = FALSE)
#' cmp <- compareEngines(ose, ~ condition,
#'                       contrast = c("condition", "tumour", "normal"))
#' cmp$overlap
#'
#' @export
compareEngines <- function(object, design, contrast = NULL, coef = NULL,
                           engines = c("DESeq2", "edgeR", "limma"),
                           alpha = 0.05, assay = NULL) {
    tables <- list()
    for (e in engines) {
        tmp <- diffExpress(object, design, contrast = contrast, coef = coef,
                           engine = e, assay = assay, name = "tmp",
                           alpha = alpha, verbose = FALSE)
        tables[[e]] <- getAnalysis(tmp, "tmp")
    }
    features <- Reduce(union, lapply(tables, `[[`, "feature"))
    merged <- data.frame(feature = features, stringsAsFactors = FALSE)
    sig <- list()
    for (e in engines) {
        tb <- tables[[e]]
        idx <- match(features, tb$feature)
        merged[[paste0("log2FC_", e)]] <- tb$log2FC[idx]
        merged[[paste0("padj_", e)]] <- tb$padj[idx]
        sig[[e]] <- tb$feature[which(tb$padj < alpha)]
    }
    consensus <- Reduce(intersect, sig)
    overlap <- c(vapply(sig, length, integer(1)),
                 consensus = length(consensus))
    list(merged = merged, overlap = overlap, consensus = consensus,
         alpha = alpha)
}
