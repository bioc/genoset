
##' @exportClass RangedDataOrGenomicRanges
setClassUnion("RangedDataOrGenomicRanges", c("RangedData", "GenomicRanges"))
##' @exportClass RangedDataOrGenoset
setClassUnion("RangedDataOrGenoset", c("RangedData", "GenoSet"))
##' @exportClass RangedDataOrGenoSetOrGenomicRanges
setClassUnion("RangedDataOrGenoSetOrGenomicRanges", c("RangedData", "GenoSet", "GenomicRanges"))



##' Get universe annotations
##'
##' Get universe annotations
##' @param x GenoSet
##' @return scalar character
##' @rdname universe-methods
##' @exportMethod universe
setMethod("universe", "GenoSet", function(x) {
  .Defunct(new="genome", msg="RangedData is being replaced with GenomicRanges. Please use the genome method.")
} )

##' @rdname universe-methods
setMethod("universe", "GRanges", function(x) {
  .Defunct(new="genome", msg="RangedData is being replaced with GenomicRanges. Please use the genome method.")
} )

##' @rdname universe-methods
##' @exportMethod "universe<-"
##' @param value scalar character, new value of universe
setMethod("universe<-", signature(x="GenoSet"),
                 function(x,value) {
                   .Defunct(new="genome", msg="RangedData is being replaced with GenomicRanges. Please use the genome method.")
                   })

##' @rdname universe-methods
setMethod("universe<-", signature(x="GRanges"),
          function(x,value) {
            .Defunct(new="genome", msg="RangedData is being replaced with GenomicRanges. Please use the genome method.")
          })

##' @rdname rownames-methods
setMethod("featureNames", signature(object="RangedData"),
          function(object) {
            .Defunct(new="rownames", msg="Please use rownames. We are switching away from eSet-specific methods.")
          })

##' @rdname rownames-methods
setMethod("featureNames<-",
          signature=signature(object="RangedData", value="ANY"),
          function(object, value) {
            .Defunct(new="rownames<-", msg="Please use rownames. We are switching away from eSet-specific methods.")
          })

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



##' Make a RangedData from segments
##'
##' Starting from a data.frame of segments, like from CBS and segTable, organize as a RangedData. Label data "score",
##' so it can easily be made into various genome browser formats using rtracklayer.
##' @param segs data.frame, like from segment in DNAcopy or segTable
##' @return RangedData
##' @family "segmented data"
##' @export 
##' @family segments
segs2RangedData <- function(segs) {
  .Defunct("segs2Granges", msg="genoset is moving towards using GenomicRanges instead of RangedData.")
}
