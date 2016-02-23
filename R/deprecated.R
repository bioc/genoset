##' @export
##' @rdname genoset-methods
setGeneric("locData<-", function(object,...,value) standardGeneric("locData<-"))
setMethod("locData<-", "GenoSet",
          function(object, ..., value) {
              .Deprecated("rowRanges<-")
              rowRanges(object) = value
              return(object)
          })

##' @export
##' @rdname genoset-methods
setGeneric("locData", function(object) standardGeneric("locData"))
setMethod("locData", "GenoSet",
          function(object) {
              .Deprecated("rowRanges")
              rowRanges(object)
          })
