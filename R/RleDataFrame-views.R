### Methods to do view operations on an RleDataFrame
###  viewMeans(Views(rledf,iranges)) works, but gives list.  Consistency good, but I really wanna vapply here ...

##' Calculate min/max/sum/mean/whichmin/whichmax over each view on each column of an RleDataFrame.
##'
##' Loop over Rles in RleDataFrame, calculate the appropriate statistic for each view. If simplify == FALSE,
##' returns a vector for each Rle. If simplify == TRUE, returns a vector for the case of a single view, otherwise,
##' a matrix. Rownames for the matrix are taken from the names of the argument \code{ir}.
##' @param x RleDataFrame
##' @param ir IRanges, views on each Rle you want to work on
##' @param na.rm scalar logical, ignore NAs?
##' @param simplify scalar logical, simplify result? For single view, vector, otherwise matrix with one row per view.
##' @param FUN S4Generic or scalar character with the name of an S4 Generic
##' @param FUN.TYPE scalar character, the storage mode for the returned vector or matrix (when simplify==TRUE).
##' @export 
##' @return With simlify == TRUE, vector for single view or single column x with simplify == TRUE, matix
##' otherwise. When simplify == FALSE, a list of length ncol(x).
##' @keywords internal
.do_rledf_views <- function(x, ir, na.rm=FALSE, simplify=TRUE, FUN, FUN.TYPE=c("numeric", "double", "integer", "logical")) {
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
setGeneric("viewSums2", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("viewSums2") })
setMethod("viewSums2", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewSums", FUN.TYPE="numeric")
          })

##' @export
setGeneric("viewMeans2", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("viewMeans2") })
setMethod("viewMeans2", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewMeans", FUN.TYPE="numeric")
          })

##' @export
setGeneric("viewMins2", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("viewMins2") })
setMethod("viewMins2", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewMins", FUN.TYPE="numeric")
          })

##' @export
setGeneric("viewMaxs2", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("viewMaxs2") })
setMethod("viewMaxs2", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewMaxs", FUN.TYPE="numeric")
          })

##' @export
setGeneric("viewWhichMins2", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("viewWhichMins2") })
setMethod("viewWhichMins2", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewWhichMins", FUN.TYPE="integer")
          })

##' @export
setGeneric("viewWhichMaxs2", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("viewWhichMaxs2") })
setMethod("viewWhichMaxs2", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, FUN="viewWhichMaxs", FUN.TYPE="integer")
          })

ds = RleDataFrame(list(a=Rle(1:5, rep(2, 5))), b=Rle(1:5, rep(2, 5)))
ir = IRanges(start=c(1, 4), end=c(3, 5), names=LETTERS[1:2])
v = Views(ds, ir)

