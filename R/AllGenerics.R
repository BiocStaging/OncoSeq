#' @title OncoSeq generics
#' @description Generic functions defining the OncoSeq accessor interface.
#'   Methods are documented with the \linkS4class{OncoSeqExperiment} class.
#' @param x An \linkS4class{OncoSeqExperiment} object.
#' @param value Replacement value.
#' @param ... Further arguments passed to methods.
#' @return Accessors return the corresponding layer; replacement functions
#'   return the modified object.
#'
#' @examples
#' ose <- exampleOncoSeq()
#' head(mutations(ose))
#' cnv(ose)
#' analyses(ose)
#' names(provenance(ose))
#'
#' @name OncoSeq-generics
NULL

#' @rdname OncoSeq-generics
#' @export
setGeneric("mutations", function(x, ...) standardGeneric("mutations"))

#' @rdname OncoSeq-generics
#' @export
setGeneric("mutations<-", function(x, value) standardGeneric("mutations<-"))

#' @rdname OncoSeq-generics
#' @export
setGeneric("cnv", function(x, ...) standardGeneric("cnv"))

#' @rdname OncoSeq-generics
#' @export
setGeneric("cnv<-", function(x, value) standardGeneric("cnv<-"))

#' @rdname OncoSeq-generics
#' @export
setGeneric("analyses", function(x, ...) standardGeneric("analyses"))

#' @rdname OncoSeq-generics
#' @export
setGeneric("analyses<-", function(x, value) standardGeneric("analyses<-"))

#' @rdname OncoSeq-generics
#' @export
setGeneric("provenance", function(x, ...) standardGeneric("provenance"))

#' @rdname OncoSeq-generics
#' @export
setGeneric("provenance<-", function(x, value) standardGeneric("provenance<-"))
