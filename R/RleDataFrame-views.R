### Methods to do view operations on an RleDataFrame
###  viewMeans(Views(rledf,iranges)) works, but gives list.  Consistency good, but I really wanna simplify via vapply here ...

### Hmm, naming is an issue. Can't add methods for viewMeans etc. without "..." in those generic function definitions.
### Maybe these shouldn't be "view" functions anyway because they aren't methods on View types.

##' @include RleDataFrame-class.R
NULL
##' Calculate min/max/sum/mean/whichmin/whichmax over each view on each column of an RleDataFrame.
##'
##' Loop over the Rle objects in an RleDataFrame, calculate the appropriate statistic for each view. If simplify == FALSE,
##' this function returns a vector for each Rle. If simplify == TRUE, it returns a vector for the case of a single view, otherwise,
##' a matrix. Rownames for the matrix are taken from the names of the argument \code{ir}.
##' @param x RleDataFrame
##' @param ir IRanges or matrix, views on every Rle. If \code{ir} is a matrix, it is converted to an IRanges using the first
##' two columns as the starts and stops. Names for the IRanges are taken from the rownames of the matrix. Such a matrix can be
##' constructed with \code{boundingIndicesByChr}.
##' @param na.rm scalar logical, ignore NAs in calculations?
##' @param simplify scalar logical, simplify result? For a single view, a vector, otherwise a matrix with one row per view.
##' @param FUN S4 Generic function or scalar character with the name of an S4 Generic function
##' @param FUN.TYPE scalar character, the storage mode for the returned vector or matrix (when simplify==TRUE).
##' @export 
##' @return With \code{simplify == TRUE}, a vector for single view or a matrix
##' otherwise. When \code{simplify == FALSE}, a list of vectors length ncol(x) where each element is of length \code{nrows(ir)}.
##' @keywords internal
##' @rdname do_rledf_views
##' @seealso RleDataFrame boundingIndicesByChr
##' @family views
.do_rledf_views <- function(x, ir, na.rm=FALSE, simplify=TRUE, FUN, FUN.TYPE=c("numeric", "double", "integer", "logical")) {
  if (is.matrix(ir)) {
    ir = IRanges(start=ir[, 1], end=ir[, 2], names=rownames(ir))
  }
  myfun = getMethod(FUN, "RleViews", where="IRanges")
  myviewfun = getMethod("Views", "Rle", where="IRanges")
  if (simplify == TRUE) {
    FUN.TYPE = match.arg(FUN.TYPE)
    nviews = length(ir)
    val = vapply(x,   	 
           FUN=function(rle) {
             myfun(myviewfun(rle, ir), na.rm=na.rm)
           }, USE.NAMES=ifelse(nviews > 1, TRUE, FALSE),
           FUN.VALUE=vector( FUN.TYPE, nviews ) )
  } else {
    val = lapply(x,   	 
      function(rle) {
        myfun(myviewfun(rle, ir), na.rm=na.rm)
      })
    }
  return(val)
}

##' @export
setGeneric("rangeSums", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeSums") })
setMethod("rangeSums", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewSums", FUN.TYPE="numeric")
          })

##' @export
setGeneric("rangeMeans", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeMeans") })
setMethod("rangeMeans", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewMeans", FUN.TYPE="numeric")
          })

##' @export
setGeneric("rangeMins", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeMins") })
setMethod("rangeMins", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewMins", FUN.TYPE="numeric")
          })

##' @export
setGeneric("rangeMaxs", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeMaxs") })
setMethod("rangeMaxs", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewMaxs", FUN.TYPE="numeric")
          })

##' @export
setGeneric("rangeWhichMins", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeWhichMins") })
setMethod("rangeWhichMins", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewWhichMins", FUN.TYPE="integer")
          })

##' @export
setGeneric("rangeWhichMaxs", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeWhichMaxs") })
setMethod("rangeWhichMaxs", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewWhichMaxs", FUN.TYPE="integer")
          })

