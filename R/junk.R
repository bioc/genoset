### All of my deprecated and defunct stuff

##' Deprecated genoset features
##'
##' GenoSet is moving towards using GenomicRanges instead of RangedData. We are also getting rid of dependencies on eSet for a potential switch to an underlying SummarizedExperiment.
##' @name genoset-deprecated
##' @aliases genoset-deprecated
NULL

##' Defunct genoset features
##'
##' The CNSet and BAFSet classes are defunct.  They only really added getter/setter methods for specific assayDataElements,
##' so they are now redundant with the preferred method of using the assayDataElement name as the third argument to bracket, e.g.
##' \code{x[i, j, "lrr"]}. Accordingly \code{BAFSet.to.ExpressionSets} is also defunct.
##'
##' Additionally, names, ranges, and space on a GenoSet are also defunct. In an effort to make a consistent API for either RangedData or
##' GRanges in the locData slot, we recommend using \code{chrNames} for \code{names} and \code{chr} for \code{space}.
##' @name genoset-defunct
##' @aliases genoset-defunct
NULL

##' @name API
##' @title Converting eSet API to RangedSummarizedExperiment
##' @param x GenoSet
##' @param object GenoSet
##' @param i name or index of assay to select
##' @param value new value to insert in assays, a matrix-like object
##' @param ... extra arguments
##' @param withDimnames scalar logical, should assays be returned with dimnames? Currently always TRUE for GenoSet.
##' @rdname api
NULL

##' @exportMethod names
setMethod("names", "GenoSet", function(x) {
  Deprecated("assayNames")
  return( assayNames(x) )
} )
