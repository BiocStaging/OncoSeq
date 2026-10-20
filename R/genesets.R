#' Read a gene-set collection from a GMT file
#'
#' @description
#' Parses a Broad-style GMT file into a named list of character vectors, the
#' format expected by [runGSEA()] and [runORA()].
#'
#' @param file Path to a `.gmt` file.
#' @param minSize,maxSize Discard sets outside this size range.
#'
#' @return A named list of gene identifiers.
#'
#' @examples
#' gmt <- system.file("extdata", "example_sets.gmt", package = "OncoSeq")
#' sets <- readGMT(gmt)
#' lengths(sets)
#'
#' @seealso [msigdbrSets()]
#' @export
readGMT <- function(file, minSize = 1L, maxSize = Inf) {
    if (!file.exists(file)) stop("file not found: ", file, call. = FALSE)
    lines <- readLines(file, warn = FALSE)
    lines <- lines[nzchar(trimws(lines))]
    parts <- strsplit(lines, "\t", fixed = TRUE)
    sets <- lapply(parts, function(p) {
        genes <- p[-c(1L, 2L)]
        unique(genes[nzchar(genes)])
    })
    names(sets) <- vapply(parts, `[`, character(1), 1L)
    sizes <- lengths(sets)
    sets[sizes >= minSize & sizes <= maxSize]
}

#' Fetch a gene-set collection through msigdbr
#'
#' @description
#' Convenience wrapper that converts an MSigDB collection obtained with the
#' `msigdbr` package into the named-list format used by OncoSeq.
#'
#' @param species Species name accepted by [msigdbr::msigdbr()], for example
#'   `"Homo sapiens"`.
#' @param collection Collection identifier, for example `"H"` for the
#'   hallmark sets or `"C2"` for curated sets.
#' @param subcollection Optional sub-collection, for example `"CP:REACTOME"`.
#'
#' @return A named list of gene symbols.
#'
#' @section Requirements:
#' Needs the `msigdbr` package, which is not a hard dependency of OncoSeq.
#'
#' @examples
#' \dontrun{
#' hallmark <- msigdbrSets("Homo sapiens", "H")
#' }
#'
#' @export
msigdbrSets <- function(species = "Homo sapiens", collection = "H",
                        subcollection = NULL) {
    .requireNamespaces("msigdbr", "msigdbrSets()")
    args <- list(species = species, collection = collection)
    if (!is.null(subcollection)) args$subcollection <- subcollection
    df <- do.call(msigdbr::msigdbr, args)
    df <- as.data.frame(df)
    nameCol <- intersect(c("gs_name", "gs_id"), colnames(df))[1L]
    geneCol <- intersect(c("gene_symbol", "ensembl_gene"), colnames(df))[1L]
    split(df[[geneCol]], df[[nameCol]])
}
