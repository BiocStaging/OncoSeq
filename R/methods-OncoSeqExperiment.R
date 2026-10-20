#' Methods for OncoSeqExperiment objects
#'
#' @description Accessors, replacement methods, subsetting and display for
#'   \linkS4class{OncoSeqExperiment}.
#'
#' @param x,object An \linkS4class{OncoSeqExperiment}.
#' @param value Replacement value.
#' @param i,j Row (feature) and column (sample) subscripts.
#' @param drop Ignored; retained for signature compatibility.
#' @param ... Further arguments, currently unused.
#'
#' @return `mutations()` a \link[S4Vectors:DataFrame-class]{DataFrame}; `cnv()` a
#'   \link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment} or `NULL`; `analyses()` a
#'   \link[S4Vectors:SimpleList-class]{SimpleList}; `provenance()` a list; `[` an
#'   \linkS4class{OncoSeqExperiment}; `show()` is called for its side effect.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' head(mutations(ose))
#' cnv(ose)
#' ose[1:5, 1:4]
#'
#' @name OncoSeqExperiment-methods
NULL

#' @rdname OncoSeqExperiment-methods
#' @export
setMethod("mutations", "OncoSeqExperiment", function(x, ...) x@mutations)

#' @rdname OncoSeqExperiment-methods
#' @export
setReplaceMethod("mutations", "OncoSeqExperiment", function(x, value) {
    x@mutations <- .normaliseMutationTable(value)
    validObject(x)
    x
})

#' @rdname OncoSeqExperiment-methods
#' @export
setMethod("cnv", "OncoSeqExperiment", function(x, ...) x@cnv)

#' @rdname OncoSeqExperiment-methods
#' @export
setReplaceMethod("cnv", "OncoSeqExperiment", function(x, value) {
    x@cnv <- .asCNVSummarizedExperiment(value)
    validObject(x)
    x
})

#' @rdname OncoSeqExperiment-methods
#' @export
setMethod("analyses", "OncoSeqExperiment", function(x, ...) x@analyses)

#' @rdname OncoSeqExperiment-methods
#' @export
setReplaceMethod("analyses", "OncoSeqExperiment", function(x, value) {
    x@analyses <- if (is(value, "SimpleList")) value else SimpleList(value)
    validObject(x)
    x
})

#' @rdname OncoSeqExperiment-methods
#' @export
setMethod("provenance", "OncoSeqExperiment", function(x, ...) x@provenance)

#' @rdname OncoSeqExperiment-methods
#' @export
setReplaceMethod("provenance", "OncoSeqExperiment", function(x, value) {
    stopifnot(is.list(value))
    x@provenance <- value
    x
})

#' @rdname OncoSeqExperiment-methods
#' @export
setMethod("[", c("OncoSeqExperiment", "ANY", "ANY", "ANY"),
    function(x, i, j, ..., drop = FALSE) {
        out <- callNextMethod(x, i, j, ..., drop = FALSE)
        if (!missing(j)) {
            keep <- colnames(out)
            if (nrow(x@mutations) > 0L) {
                out@mutations <- x@mutations[
                    as.character(x@mutations$sample) %in% keep, , drop = FALSE]
            }
            if (!is.null(x@cnv)) {
                out@cnv <- x@cnv[, intersect(colnames(x@cnv), keep),
                                 drop = FALSE]
            }
        }
        out
    })

#' @rdname OncoSeqExperiment-methods
#' @importFrom methods show
#' @export
setMethod("show", "OncoSeqExperiment", function(object) {
    cat("class: OncoSeqExperiment\n")
    cat("dim:", nrow(object), "features x", ncol(object), "samples\n")
    cat("assays(", length(SummarizedExperiment::assayNames(object)), "): ",
        paste(SummarizedExperiment::assayNames(object), collapse = ", "),
        "\n", sep = "")
    cat(.showLine("colData", colnames(SummarizedExperiment::colData(object))))
    cat(.showLine("rowData", colnames(SummarizedExperiment::rowData(object))))

    mut <- object@mutations
    if (nrow(mut) > 0L) {
        cat("mutations: ", nrow(mut), " variants in ",
            length(unique(as.character(mut$gene))), " genes across ",
            length(unique(as.character(mut$sample))), " samples\n", sep = "")
    } else {
        cat("mutations: none\n")
    }

    if (!is.null(object@cnv)) {
        cat("cnv: ", nrow(object@cnv), " genes x ", ncol(object@cnv),
            " samples\n", sep = "")
    } else {
        cat("cnv: none\n")
    }

    if (length(object@analyses) > 0L) {
        cat("analyses(", length(object@analyses), "): ",
            paste(names(object@analyses), collapse = ", "), "\n", sep = "")
    } else {
        cat("analyses: none\n")
    }
})
