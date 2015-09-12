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

##' @exportMethod pData
##' @rdname api
setMethod("pData",
          signature=signature(object="GenoSet"),
          function(object) {
              .Deprecated("colData")
              object@phenoData@data
          })

##' @exportMethod colData
##' @rdname api
setMethod("colData",
          signature=signature(x="GenoSet"),
          function(x) {
              x@phenoData@data
          })

##' @rdname rownames-methods
##' @exportMethod featureNames
setMethod("featureNames", signature(object="GenomicRanges"),
          function(object) {
              .Deprecated("rownames")
            names(object)
          })

##' @rdname rownames-methods
##' @param object GenoSet
##' @exportMethod featureNames
setMethod("featureNames", signature(object="GenoSet"),
          function(object) {
              .Deprecated("rownames")
            rownames(object)
          })

##' @exportMethod "featureNames<-"
##' @rdname rownames-methods
setMethod("featureNames<-",
          signature=signature(object="GenomicRanges", value="ANY"),
          function(object, value) {
              .Deprecated("rownames<-")
              names(object) = value
            return(object)
          })

##' @rdname locData-methods
##' @exportMethod rowRanges
##' @importFrom SummarizedExperiment rowRanges
##' @param x GenoSet object
setMethod("rowRanges",
          signature=signature(x="GenoSet"),
          function(x) {
              locData(x)
          })

##' @rdname colnames
##' @param object a Genoset
##' @exportMethod sampleNames
setMethod("sampleNames", signature(object="GenoSet"),
          function(object) {
              .Deprecated("colnames")
            colnames(object)
          })

##' @exportMethod "sampleNames<-"
##' @rdname colnames
setMethod("sampleNames<-", signature(object="GenoSet"),
          function(object, value) {
              .Deprecated("colnames<-")
            colnames(object) = value
            return(object)
          })

##' @exportMethod assayData
##' @rdname api
setMethod("assays", signature(x="GenoSet"),
          function(x, ..., withDimnames=TRUE) {
              x@assayData
          })

##' @exportMethod assay
##' @rdname api
setMethod("assay", signature(x="GenoSet",i="ANY"),
          function(x,i) {
              x@assayData[[i]]
          })


##' @exportMethod "assay<-"
##' @rdname api
setMethod("assay<-", signature(x="GenoSet",i="ANY",value="ANY"),
          function(x,i,value) {
              assayDataElement(x,i) <- value
              return(x)
          })

##' @exportMethod assayNames
##' @rdname api
setMethod("assayNames", signature(x="GenoSet"),
          function(x) {
              names(x@assayData)
          })
