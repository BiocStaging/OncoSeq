#' @importFrom survival Surv survfit survdiff coxph
#' @importFrom stats pchisq setNames
NULL

#' Kaplan-Meier analysis of a gene or clinical variable
#'
#' @description
#' Splits samples into groups by expression, mutation status, copy-number
#' state or a clinical variable, fits Kaplan-Meier curves and tests them
#' with the log-rank statistic, and reports the hazard ratio from the
#' corresponding Cox model.
#'
#' @details
#' For expression the split may be at the median, by tertiles, by quartiles
#' (comparing the outer two), or at an optimal cut-point. The optimal
#' cut-point scans candidate thresholds between the 25th and 75th
#' percentiles and keeps the one maximising the log-rank statistic. That
#' search inflates significance, so the reported `p_logrank_adjusted` applies
#' the Altman approximation; the unadjusted value is reported alongside it
#' and should not be quoted on its own.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param feature Gene identifier or `colData()` column name to stratify on.
#' @param by Layer the feature belongs to: `"expression"`, `"mutation"`,
#'   `"cnv"`, `"clinical"`, or `"auto"` to search them in that order.
#' @param split Grouping rule for continuous features: `"median"`,
#'   `"tertile"`, `"quartile"` or `"optimal"`.
#' @param timeColumn,eventColumn `colData()` columns holding follow-up time
#'   and the event indicator (1 = event, 0 = censored).
#' @param assay Expression assay to use; log-scale values are derived from
#'   counts when needed.
#' @param minGroupSize Refuse to test when a group is smaller than this.
#'
#' @return A list of class `oncoseq_survival` with the `survfit` object, the
#'   grouping factor, the log-rank test, the Cox hazard ratio with
#'   confidence interval, and the cut-point when one was chosen.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' fit <- survivalFit(ose, "TP53", by = "mutation")
#' fit$p_logrank
#'
#' @seealso [coxModel()], [survivalScreen()], [plotKM()]
#' @export
survivalFit <- function(object, feature,
                        by = c("auto", "expression", "mutation", "cnv",
                               "clinical"),
                        split = c("median", "tertile", "quartile",
                                  "optimal"),
                        timeColumn = "os_time", eventColumn = "os_event",
                        assay = NULL, minGroupSize = 5L) {
    by <- match.arg(by)
    split <- match.arg(split)
    surv <- .survivalResponse(object, timeColumn, eventColumn)

    grouping <- .featureGrouping(object, feature, by, split, assay, surv)
    group <- grouping$group
    df <- data.frame(time = surv$time, event = surv$event, group = group)
    df <- df[!is.na(df$group), , drop = FALSE]
    df$group <- droplevels(factor(df$group))

    sizes <- table(df$group)
    if (length(sizes) < 2L) {
        stop("stratifying on ", sQuote(feature), " produced a single group",
             call. = FALSE)
    }
    if (any(sizes < minGroupSize)) {
        stop("group(s) ", paste(names(sizes)[sizes < minGroupSize],
             collapse = ", "), " have fewer than ", minGroupSize,
             " samples; use a different split", call. = FALSE)
    }

    fit <- survival::survfit(survival::Surv(time, event) ~ group, data = df)
    sd <- survival::survdiff(survival::Surv(time, event) ~ group, data = df)
    pLogrank <- stats::pchisq(sd$chisq, df = length(sizes) - 1L,
                              lower.tail = FALSE)
    cox <- survival::coxph(survival::Surv(time, event) ~ group, data = df)
    cs <- summary(cox)

    out <- list(
        fit = fit,
        data = df,
        feature = feature,
        layer = grouping$layer,
        split = if (grouping$continuous) split else NA_character_,
        cutpoint = grouping$cutpoint,
        n = as.integer(sizes),
        groups = names(sizes),
        chisq = unname(sd$chisq),
        p_logrank = pLogrank,
        p_logrank_adjusted = if (identical(split, "optimal") &&
                                 grouping$continuous) {
            .altmanAdjust(pLogrank)
        } else pLogrank,
        hazard_ratio = unname(cs$conf.int[, "exp(coef)"]),
        hr_lower = unname(cs$conf.int[, "lower .95"]),
        hr_upper = unname(cs$conf.int[, "upper .95"]),
        cox = cox)
    class(out) <- "oncoseq_survival"
    out
}

#' @export
print.oncoseq_survival <- function(x, ...) {
    cat("Kaplan-Meier analysis of", x$feature,
        paste0("(", x$layer, ")"), "\n")
    cat("groups:", paste(paste0(x$groups, " (n=", x$n, ")"),
                         collapse = ", "), "\n")
    if (!is.na(x$cutpoint)) cat("cut-point:", signif(x$cutpoint, 4), "\n")
    cat("log-rank p:", format.pval(x$p_logrank, digits = 3), "\n")
    if (!isTRUE(all.equal(x$p_logrank, x$p_logrank_adjusted))) {
        cat("log-rank p, adjusted for cut-point search:",
            format.pval(x$p_logrank_adjusted, digits = 3), "\n")
    }
    cat("hazard ratio:", paste(signif(x$hazard_ratio, 3), collapse = ", "),
        "\n")
    invisible(x)
}

## Altman's approximation for a cut-point chosen by maximising the log-rank
## statistic over the inner 50% of the distribution.
.altmanAdjust <- function(p) {
    if (is.na(p) || p <= 0) return(p)
    z <- stats::qnorm(1 - p / 2)
    eps <- 0.25
    adj <- stats::dnorm(z) * (z - 1 / z) *
        log(((1 - eps)^2) / (eps^2)) + 4 * stats::dnorm(z) / z
    min(1, max(p, adj))
}

.survivalResponse <- function(object, timeColumn, eventColumn) {
    cd <- SummarizedExperiment::colData(object)
    for (nm in c(timeColumn, eventColumn)) {
        if (!nm %in% colnames(cd)) {
            stop("colData has no column ", sQuote(nm),
                 ". Survival analysis needs follow-up time and an event ",
                 "indicator; available columns: ",
                 paste(colnames(cd), collapse = ", "), call. = FALSE)
        }
    }
    time <- suppressWarnings(as.numeric(cd[[timeColumn]]))
    event <- cd[[eventColumn]]
    if (is.logical(event)) event <- as.integer(event)
    if (is.factor(event) || is.character(event)) {
        event <- suppressWarnings(as.integer(sub(":.*$", "",
                                                 as.character(event))))
    }
    event <- as.integer(event)
    bad <- !event %in% c(0L, 1L, NA_integer_)
    if (any(bad, na.rm = TRUE)) {
        stop(sQuote(eventColumn), " must be coded 0 = censored, 1 = event",
             call. = FALSE)
    }
    keep <- !is.na(time) & !is.na(event) & time >= 0
    if (sum(keep) < 10L) {
        stop("only ", sum(keep), " sample(s) have complete follow-up data",
             call. = FALSE)
    }
    time[!keep] <- NA_real_
    event[!keep] <- NA_integer_
    list(time = time, event = event, keep = keep)
}

.featureGrouping <- function(object, feature, by, split, assay, surv) {
    cd <- SummarizedExperiment::colData(object)
    tryOrder <- if (by == "auto") {
        c("expression", "mutation", "cnv", "clinical")
    } else by

    for (layer in tryOrder) {
        if (layer == "expression") {
            mat <- .logExpressionMatrix(object, assay)
            if (feature %in% rownames(mat)) {
                v <- as.numeric(mat[feature, ])
                g <- .splitContinuous(v, split, surv)
                return(list(group = g$group, layer = "expression",
                            continuous = TRUE, cutpoint = g$cutpoint))
            }
        } else if (layer == "mutation") {
            mut <- mutations(object)
            if (nrow(mut) > 0L && feature %in% as.character(mut$gene)) {
                hit <- unique(as.character(
                    mut$sample[as.character(mut$gene) == feature]))
                g <- factor(ifelse(colnames(object) %in% hit,
                                   "mutant", "wild-type"),
                            levels = c("wild-type", "mutant"))
                return(list(group = g, layer = "mutation",
                            continuous = FALSE, cutpoint = NA_real_))
            }
        } else if (layer == "cnv") {
            cn <- cnv(object)
            if (!is.null(cn) && feature %in% rownames(cn)) {
                v <- as.numeric(SummarizedExperiment::assay(cn, 1L)[feature, ])
                names(v) <- colnames(cn)
                v <- v[match(colnames(object), names(v))]
                g <- factor(ifelse(is.na(v), NA_character_,
                            ifelse(v >= 1, "gain",
                            ifelse(v <= -1, "loss", "neutral"))),
                            levels = c("loss", "neutral", "gain"))
                return(list(group = droplevels(g), layer = "cnv",
                            continuous = FALSE, cutpoint = NA_real_))
            }
        } else if (layer == "clinical") {
            if (feature %in% colnames(cd)) {
                v <- cd[[feature]]
                if (is.numeric(v)) {
                    g <- .splitContinuous(as.numeric(v), split, surv)
                    return(list(group = g$group, layer = "clinical",
                                continuous = TRUE, cutpoint = g$cutpoint))
                }
                return(list(group = factor(as.character(v)),
                            layer = "clinical", continuous = FALSE,
                            cutpoint = NA_real_))
            }
        }
    }
    where <- if (by == "auto") {
        "expression, mutation, copy-number or clinical data"
    } else {
        paste0(by, " layer")
    }
    stop("feature ", sQuote(feature), " was not found in the ", where,
         call. = FALSE)
}

.splitContinuous <- function(v, split, surv) {
    vk <- v
    vk[!surv$keep] <- NA_real_
    cutpoint <- NA_real_
    if (split == "median") {
        cutpoint <- stats::median(vk, na.rm = TRUE)
        group <- factor(ifelse(vk > cutpoint, "high", "low"),
                        levels = c("low", "high"))
    } else if (split == "tertile") {
        qs <- stats::quantile(vk, c(1/3, 2/3), na.rm = TRUE)
        group <- cut(vk, breaks = c(-Inf, qs, Inf),
                     labels = c("low", "mid", "high"))
    } else if (split == "quartile") {
        qs <- stats::quantile(vk, c(0.25, 0.75), na.rm = TRUE)
        group <- factor(ifelse(vk <= qs[1L], "low",
                        ifelse(vk >= qs[2L], "high", NA_character_)),
                        levels = c("low", "high"))
    } else {
        cutpoint <- .optimalCutpoint(vk, surv$time, surv$event)
        group <- factor(ifelse(vk > cutpoint, "high", "low"),
                        levels = c("low", "high"))
    }
    list(group = group, cutpoint = cutpoint, values = vk)
}

## Scan candidate thresholds over the inner half of the distribution and
## keep the one giving the largest log-rank chi-square. The search inflates
## significance; survivalFit() reports an Altman-corrected p-value with it.
.optimalCutpoint <- function(v, time, event, minFraction = 0.25) {
    ok <- !is.na(v) & !is.na(time) & !is.na(event)
    vv <- v[ok]
    tt <- time[ok]
    ee <- event[ok]
    bounds <- stats::quantile(vv, c(minFraction, 1 - minFraction),
                              na.rm = TRUE)
    candidates <- sort(unique(vv[vv >= bounds[1L] & vv <= bounds[2L]]))
    candidates <- candidates[-length(candidates)]
    if (length(candidates) == 0L) {
        return(stats::median(vv, na.rm = TRUE))
    }
    chisq <- vapply(candidates, function(cp) {
        g <- vv > cp
        if (length(unique(g)) < 2L) return(NA_real_)
        sd <- try(survival::survdiff(survival::Surv(tt, ee) ~ g),
                  silent = TRUE)
        if (inherits(sd, "try-error")) NA_real_ else unname(sd$chisq)
    }, numeric(1))
    if (all(is.na(chisq))) return(stats::median(vv, na.rm = TRUE))
    candidates[which.max(chisq)]
}

#' Multivariable Cox proportional-hazards model
#'
#' @description
#' Fits [survival::coxph()] with covariates drawn from `colData()`, from
#' expression (prefix a gene name with `expr:`) or from mutation status
#' (prefix with `mut:`), and returns a tidy table of hazard ratios.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param covariates Character vector of terms, for example
#'   `c("age", "stage", "expr:TP53", "mut:RB1")`.
#' @param timeColumn,eventColumn Follow-up and event columns in `colData()`.
#' @param assay Expression assay used for `expr:` terms.
#' @param name Name under which the result is stored in the registry; `NULL`
#'   returns the table instead of the object.
#'
#' @return The object with the tidy table stored, or the table itself when
#'   `name` is `NULL`. Columns are `term`, `coef`, `HR`, `ci_lower`,
#'   `ci_upper`, `se`, `z` and `pvalue`. The fitted model is attached as the
#'   `model` attribute.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' coxModel(ose, c("expr:TP53", "mut:TP53"), name = NULL)
#'
#' @export
coxModel <- function(object, covariates, timeColumn = "os_time",
                     eventColumn = "os_event", assay = NULL,
                     name = "cox") {
    surv <- .survivalResponse(object, timeColumn, eventColumn)
    df <- data.frame(time = surv$time, event = surv$event)
    terms <- character(0)
    for (cv in covariates) {
        res <- .covariateVector(object, cv, assay)
        df[[res$name]] <- res$value
        terms <- c(terms, res$name)
    }
    df <- df[stats::complete.cases(df), , drop = FALSE]
    if (nrow(df) < 10L) {
        stop("only ", nrow(df), " samples have complete data for all ",
             "covariates", call. = FALSE)
    }
    fml <- stats::as.formula(paste("survival::Surv(time, event) ~",
                                   paste(terms, collapse = " + ")))
    fit <- survival::coxph(fml, data = df)
    s <- summary(fit)
    out <- data.frame(
        term = rownames(s$coefficients),
        coef = s$coefficients[, "coef"],
        HR = s$conf.int[, "exp(coef)"],
        ci_lower = s$conf.int[, "lower .95"],
        ci_upper = s$conf.int[, "upper .95"],
        se = s$coefficients[, "se(coef)"],
        z = s$coefficients[, "z"],
        pvalue = s$coefficients[, "Pr(>|z|)"],
        stringsAsFactors = FALSE)
    rownames(out) <- NULL
    attr(out, "model") <- fit
    attr(out, "n") <- nrow(df)
    attr(out, "n_events") <- sum(df$event)
    if (is.null(name)) return(out)
    object <- .addAnalysis(object, name, out, analysis = "Cox regression",
        covariates = covariates, n_samples = nrow(df))
    .addProvenance(object, "coxModel", name = name)
}

.covariateVector <- function(object, spec, assay) {
    cd <- SummarizedExperiment::colData(object)
    if (grepl("^expr:", spec)) {
        gene <- sub("^expr:", "", spec)
        mat <- .logExpressionMatrix(object, assay)
        if (!gene %in% rownames(mat)) {
            stop("gene ", sQuote(gene), " is not in the expression data",
                 call. = FALSE)
        }
        return(list(name = make.names(spec), value = as.numeric(mat[gene, ])))
    }
    if (grepl("^mut:", spec)) {
        gene <- sub("^mut:", "", spec)
        mut <- mutations(object)
        hit <- unique(as.character(mut$sample[as.character(mut$gene) == gene]))
        return(list(name = make.names(spec),
                    value = factor(ifelse(colnames(object) %in% hit,
                                          "mutant", "wild-type"),
                                   levels = c("wild-type", "mutant"))))
    }
    if (!spec %in% colnames(cd)) {
        stop("covariate ", sQuote(spec), " is not a colData column; prefix ",
             "gene names with 'expr:' or 'mut:'", call. = FALSE)
    }
    v <- cd[[spec]]
    if (is.character(v)) v <- factor(v)
    list(name = make.names(spec), value = v)
}

#' Genome-wide screen for survival-associated genes
#'
#' @description
#' Fits a Cox proportional-hazards model for every gene, with expression as
#' a continuous covariate and optional clinical adjustment, and corrects the
#' resulting p-values across genes with Benjamini-Hochberg.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param genes Genes to test; defaults to all features.
#' @param covariates Optional clinical covariates to adjust for.
#' @param timeColumn,eventColumn Follow-up and event columns.
#' @param assay Expression assay; log-scale values are derived from counts
#'   when needed.
#' @param scale Standardise each gene before fitting, so that the hazard
#'   ratio refers to a one standard-deviation increase and is comparable
#'   across genes.
#' @param minVariance Skip genes whose variance is at or below this value.
#' @param name Name under which the result is stored; `NULL` returns the
#'   table.
#'
#' @return The object with a table of `gene`, `HR`, `ci_lower`, `ci_upper`,
#'   `z`, `pvalue` and `padj`, ordered by adjusted p-value.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' res <- survivalScreen(ose, genes = rownames(ose)[1:50], name = NULL)
#' head(res)
#'
#' @export
survivalScreen <- function(object, genes = NULL, covariates = NULL,
                           timeColumn = "os_time", eventColumn = "os_event",
                           assay = NULL, scale = TRUE, minVariance = 0,
                           name = "survival_screen") {
    surv <- .survivalResponse(object, timeColumn, eventColumn)
    mat <- .logExpressionMatrix(object, assay)
    if (!is.null(genes)) {
        missing_genes <- setdiff(genes, rownames(mat))
        if (length(missing_genes) == length(genes)) {
            stop("none of the requested genes is in the expression data",
                 call. = FALSE)
        }
        mat <- mat[intersect(genes, rownames(mat)), , drop = FALSE]
    }

    base <- data.frame(time = surv$time, event = surv$event)
    extraTerms <- character(0)
    for (cv in covariates %||% character(0)) {
        res <- .covariateVector(object, cv, assay)
        base[[res$name]] <- res$value
        extraTerms <- c(extraTerms, res$name)
    }

    rhs <- paste(c("gene", extraTerms), collapse = " + ")
    fml <- stats::as.formula(paste("survival::Surv(time, event) ~", rhs))

    vars <- matrixStats::rowVars(mat, na.rm = TRUE)
    testable <- which(vars > minVariance)

    rows <- lapply(testable, function(i) {
        g <- as.numeric(mat[i, ])
        if (scale) g <- as.numeric(scale(g))
        d <- base
        d$gene <- g
        d <- d[stats::complete.cases(d), , drop = FALSE]
        if (nrow(d) < 10L || sum(d$event) < 3L) return(NULL)
        fit <- try(survival::coxph(fml, data = d), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        s <- summary(fit)
        if (!"gene" %in% rownames(s$coefficients)) return(NULL)
        data.frame(
            gene = rownames(mat)[i],
            HR = s$conf.int["gene", "exp(coef)"],
            ci_lower = s$conf.int["gene", "lower .95"],
            ci_upper = s$conf.int["gene", "upper .95"],
            z = s$coefficients["gene", "z"],
            pvalue = s$coefficients["gene", "Pr(>|z|)"],
            n = nrow(d),
            stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
    if (is.null(out) || nrow(out) == 0L) {
        stop("no gene could be fitted; check follow-up data and gene ",
             "variance", call. = FALSE)
    }
    out$padj <- .bh(out$pvalue)
    out <- out[order(out$padj, out$pvalue, na.last = TRUE), , drop = FALSE]
    rownames(out) <- NULL

    if (is.null(name)) return(out)
    object <- .addAnalysis(object, name, out, analysis = "survival screen",
        covariates = covariates, nTested = nrow(out))
    .addProvenance(object, "survivalScreen", name = name)
}
