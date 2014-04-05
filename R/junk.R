
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
##' @exportMethod featureNames
setMethod("featureNames", signature(object="GenoSet"),
          function(object) {
             .Defunct(new="rownames", msg="Please use rownames. We are switching away from eSet-specific methods.")
          })

##' @rdname rownames-methods
setMethod("featureNames", signature(object="GRanges"),
          function(object) {
            .Defunct(new="rownames", msg="Please use rownames. We are switching away from eSet-specific methods.")
          })

##' @rdname rownames-methods
setMethod("featureNames<-",
          signature=signature(object="GenoSet", value="ANY"),
          function(object, value) {
            .Defunct(new="rownames<-", msg="Please use rownames. We are switching away from eSet-specific methods.")
          })
##' @rdname rownames-methods
setMethod("featureNames<-",
          signature=signature(object="GRanges", value="ANY"),
          function(object, value) {
            .Defunct(new="rownames<-", msg="Please use rownames. We are switching away from eSet-specific methods.")
          })
##' @rdname rownames-methods
setMethod("featureNames<-",
          signature=signature(object="RangedData", value="ANY"),
          function(object, value) {
            .Defunct(new="rownames<-", msg="Please use rownames. We are switching away from eSet-specific methods.")
          })
