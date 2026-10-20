#' @importFrom stats binom.test
NULL

## Variant classifications that alter the protein product. Used as the
## default filter for driver-gene frequency and mutation burden, matching
## the convention of the TCGA MC3 and maftools analyses.
NONSYNONYMOUS_CLASSES <- c(
    "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation",
    "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
    "Splice_Site", "Translation_Start_Site",
    ## cBioPortal spellings
    "Missense_Mutation", "Nonsense_Mutation", "Splice_Region",
    "Frameshift_Deletion", "Frameshift_Insertion")

#' Read a MAF file
#'
#' @description
#' Reads a Mutation Annotation Format file (optionally gzip-compressed) and
#' harmonises its column names to the OncoSeq vocabulary, keeping all
#' original columns.
#'
#' @param file Path to a MAF or MAF-like delimited file.
#' @param sep Field separator.
#' @param comment Lines beginning with this character are skipped; MAF files
#'   from the GDC begin with a version header.
#'
#' @return A \link[S4Vectors:DataFrame-class]{DataFrame} with at least `sample`, `gene` and,
#'   when present, `variant_class`.
#'
#' @examples
#' maf <- system.file("extdata", "example_mutations.maf", package = "OncoSeq")
#' head(importMAF(maf))
#'
#' @seealso [mutationFrequency()], [oncoplotData()]
#' @export
importMAF <- function(file, sep = "\t", comment = "#") {
    .normaliseMutationTable(.readMAFTable(file, sep = sep, comment = comment))
}

.readMAFTable <- function(file, sep = "\t", comment = "#") {
    if (!file.exists(file)) stop("file not found: ", file, call. = FALSE)
    df <- utils::read.delim(file, sep = sep, header = TRUE,
                            comment.char = comment, quote = "",
                            check.names = FALSE, stringsAsFactors = FALSE)
    if (ncol(df) < 2L) {
        stop("only one column was parsed from ", file,
             "; check the 'sep' argument", call. = FALSE)
    }
    df
}

#' Summarise the somatic mutation layer
#'
#' @description
#' Produces per-sample and per-gene summaries of the mutation table: the
#' number of variants, the number of distinct genes or samples affected, and
#' the breakdown by variant classification.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param nonsynonymousOnly Restrict to protein-altering variants.
#'
#' @return A list with elements `perSample`, `perGene` and `byClass`, each a
#'   `data.frame`.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' s <- mutationSummary(ose)
#' head(s$perGene)
#'
#' @export
mutationSummary <- function(object, nonsynonymousOnly = TRUE) {
    mut <- .mutationTable(object, nonsynonymousOnly)
    samples <- colnames(object)

    perSample <- data.frame(
        sample = samples,
        n_variants = as.integer(table(factor(mut$sample, levels = samples))),
        stringsAsFactors = FALSE)
    perSample$n_genes <- vapply(samples, function(s) {
        length(unique(mut$gene[mut$sample == s]))
    }, integer(1))

    genes <- sort(unique(mut$gene))
    perGene <- data.frame(
        gene = genes,
        n_variants = as.integer(table(factor(mut$gene, levels = genes))),
        stringsAsFactors = FALSE)
    perGene$n_samples <- vapply(genes, function(g) {
        length(unique(mut$sample[mut$gene == g]))
    }, integer(1))
    perGene <- perGene[order(-perGene$n_samples, -perGene$n_variants), ,
                       drop = FALSE]
    rownames(perGene) <- NULL

    byClass <- if ("variant_class" %in% colnames(mut)) {
        tab <- table(mut$variant_class)
        data.frame(variant_class = names(tab),
                   n = as.integer(tab),
                   stringsAsFactors = FALSE)
    } else {
        data.frame(variant_class = character(0), n = integer(0))
    }
    byClass <- byClass[order(-byClass$n), , drop = FALSE]
    rownames(byClass) <- NULL

    list(perSample = perSample, perGene = perGene, byClass = byClass)
}

#' Mutation frequency across the cohort
#'
#' @description
#' For each gene, the proportion of samples carrying at least one mutation,
#' with an exact binomial (Clopper-Pearson) confidence interval. The
#' denominator is the number of samples in the object, so genes absent from
#' the mutation table are correctly treated as unmutated.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param genes Optional character vector restricting the output.
#' @param topN Return only the `topN` most frequently mutated genes; `Inf`
#'   returns all.
#' @param nonsynonymousOnly Restrict to protein-altering variants.
#' @param conf.level Confidence level for the interval.
#'
#' @return A `data.frame` with `gene`, `n_mutated`, `n_samples`,
#'   `frequency`, `ci_lower` and `ci_upper`, ordered by frequency.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' mutationFrequency(ose, topN = 5)
#'
#' @export
mutationFrequency <- function(object, genes = NULL, topN = 20,
                              nonsynonymousOnly = TRUE, conf.level = 0.95) {
    mut <- .mutationTable(object, nonsynonymousOnly)
    nSamples <- ncol(object)
    if (!is.null(genes)) mut <- mut[mut$gene %in% genes, , drop = FALSE]

    if (nrow(mut) == 0L) {
        return(data.frame(gene = character(0), n_mutated = integer(0),
                          n_samples = integer(0), frequency = numeric(0),
                          ci_lower = numeric(0), ci_upper = numeric(0)))
    }

    counts <- tapply(mut$sample, mut$gene,
                     function(x) length(unique(x)))
    df <- data.frame(
        gene = names(counts),
        n_mutated = as.integer(counts),
        n_samples = nSamples,
        stringsAsFactors = FALSE)
    df$frequency <- df$n_mutated / nSamples

    ci <- vapply(df$n_mutated, function(k) {
        stats::binom.test(k, nSamples, conf.level = conf.level)$conf.int
    }, numeric(2))
    df$ci_lower <- ci[1L, ]
    df$ci_upper <- ci[2L, ]

    df <- df[order(-df$frequency, df$gene), , drop = FALSE]
    rownames(df) <- NULL
    if (is.finite(topN)) df <- utils::head(df, topN)
    df
}

#' Tumour mutation burden
#'
#' @description
#' Non-synonymous somatic variants per megabase of sequenced territory for
#' every sample. The default capture size of 38 Mb is the figure
#' conventionally used for whole-exome panels; supply the actual footprint
#' of your assay when it is known.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param captureSize Size of the sequenced territory in megabases.
#' @param nonsynonymousOnly Restrict to protein-altering variants.
#' @param store Add a `tmb` column to `colData()` and return the object.
#'
#' @return The object with a `tmb` column, or a `data.frame` of burdens when
#'   `store = FALSE`.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' tmb <- tumourMutationBurden(ose, store = FALSE)
#' summary(tmb$tmb)
#'
#' @export
tumourMutationBurden <- function(object, captureSize = 38,
                                 nonsynonymousOnly = TRUE, store = TRUE) {
    if (captureSize <= 0) {
        stop("'captureSize' must be a positive number of megabases",
             call. = FALSE)
    }
    mut <- .mutationTable(object, nonsynonymousOnly)
    samples <- colnames(object)
    n <- as.integer(table(factor(mut$sample, levels = samples)))
    df <- data.frame(sample = samples, n_nonsynonymous = n,
                     tmb = n / captureSize, stringsAsFactors = FALSE)
    if (!store) return(df)
    SummarizedExperiment::colData(object)$tmb <- df$tmb
    SummarizedExperiment::colData(object)$n_nonsynonymous <-
        df$n_nonsynonymous
    .addProvenance(object, "tumourMutationBurden", captureSize = captureSize)
}

#' Build the matrix behind an oncoprint
#'
#' @description
#' Returns a gene-by-sample character matrix whose entries are the variant
#' classification carried by that gene in that sample, empty where the gene
#' is not mutated. Samples carrying several variants in one gene are
#' labelled `"Multi_Hit"`, the convention used by oncoprint renderers.
#'
#' @param object An \linkS4class{OncoSeqExperiment}.
#' @param genes Genes to include; defaults to the `topN` most frequently
#'   mutated.
#' @param topN Number of genes when `genes` is `NULL`.
#' @param nonsynonymousOnly Restrict to protein-altering variants.
#' @param sortSamples Order samples by mutation load, so the staircase
#'   pattern of an oncoprint emerges.
#'
#' @return A character matrix with genes in rows and samples in columns,
#'   ready for [ComplexHeatmap::oncoPrint()] or [plotOncoprint()].
#'
#' @examples
#' ose <- exampleOncoSeq()
#' m <- oncoplotData(ose, topN = 5)
#' m[, 1:4]
#'
#' @export
oncoplotData <- function(object, genes = NULL, topN = 20,
                         nonsynonymousOnly = TRUE, sortSamples = TRUE) {
    mut <- .mutationTable(object, nonsynonymousOnly)
    if (is.null(genes)) {
        genes <- mutationFrequency(object, topN = topN,
            nonsynonymousOnly = nonsynonymousOnly)$gene
    }
    if (length(genes) == 0L) {
        stop("no mutated genes available to plot", call. = FALSE)
    }
    samples <- colnames(object)
    mat <- matrix("", nrow = length(genes), ncol = length(samples),
                  dimnames = list(genes, samples))

    mut <- mut[mut$gene %in% genes, , drop = FALSE]
    cls <- if ("variant_class" %in% colnames(mut)) {
        as.character(mut$variant_class)
    } else rep("Mutated", nrow(mut))

    key <- paste(mut$gene, mut$sample, sep = "\r")
    first <- !duplicated(key)
    multi <- key %in% key[duplicated(key)]
    value <- ifelse(multi, "Multi_Hit", cls)

    idx <- cbind(match(mut$gene[first], genes),
                 match(mut$sample[first], samples))
    mat[idx] <- value[first]

    if (sortSamples) {
        load <- colSums(mat != "")
        mat <- mat[, order(-load), drop = FALSE]
    }
    mat
}

.mutationTable <- function(object, nonsynonymousOnly = TRUE) {
    stopifnot(is(object, "OncoSeqExperiment"))
    mut <- as.data.frame(mutations(object))
    if (nrow(mut) == 0L) {
        stop("this object carries no mutation data; attach a variant table ",
             "with mutations(object) <- ... or import one with importMAF()",
             call. = FALSE)
    }
    if (nonsynonymousOnly && "variant_class" %in% colnames(mut)) {
        keep <- mut$variant_class %in% NONSYNONYMOUS_CLASSES
        if (!any(keep)) {
            warning("no variant matches the non-synonymous classification ",
                    "list; using all variants. Observed classes: ",
                    paste(utils::head(unique(mut$variant_class), 5L),
                          collapse = ", "), call. = FALSE)
        } else {
            mut <- mut[keep, , drop = FALSE]
        }
    }
    mut
}
