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
##' @param RLEFUN function, internal rle view summary function like .rle_view_means
##' @param FUN.TYPE scalar character, the storage mode for the returned vector or matrix (when simplify==TRUE).
##' @return With \code{simplify == TRUE}, a vector for single view or a matrix
##' otherwise. When \code{simplify == FALSE}, a list of vectors length ncol(x) where each element is of length \code{nrows(ir)}.
##' @keywords internal
##' @rdname do_rledf_views
##' @seealso RleDataFrame boundingIndicesByChr
##' @family views
.do_rledf_views <- function(x, ir, na.rm=FALSE, simplify=TRUE, RLEFUN, FUN.TYPE=c("numeric", "double", "integer", "logical")) {
  # Make an IRanges from ranges matrix if necessary
  if (is.matrix(ir)) {
    ir = IRanges(start=ir[, 1], end=ir[, 2], names=rownames(ir))
  }
  # Trim IRanges once if necessary
  start(ir)[ start(ir) < 1L ] = 1L
  end(ir)[ end(ir) > nrow(x) ] = nrow(x)
  # Hoist the Views dispatch
  myviewfun = getMethod("Views", "Rle", where="IRanges")
  # Calculate the view stats
  if (simplify == TRUE) {
    FUN.TYPE = match.arg(FUN.TYPE)
    nviews = length(ir)
    val = vapply(x,   	 
           FUN=function(rle) {
             RLEFUN(myviewfun(rle, ir), na.rm=na.rm)
           }, USE.NAMES=nviews > 1, 
           FUN.VALUE=vector( FUN.TYPE, nviews ) )
  } else {
    val = lapply(x,   	 
      function(rle) {
        RLEFUN(myviewfun(rle, ir), na.rm=na.rm)
      })
    }
  return(val)
}

##' @export rangeSums
setGeneric("rangeSums", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeSums") })
setMethod("rangeSums", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_sums, FUN.TYPE="numeric")
          })

##' @export rangeMeans
setGeneric("rangeMeans", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeMeans") })
setMethod("rangeMeans", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_means, FUN.TYPE="numeric")
          })

##' @export rangeMins
setGeneric("rangeMins", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeMins") })
setMethod("rangeMins", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_mins, FUN.TYPE="numeric")
          })

##' @export rangeMaxs
setGeneric("rangeMaxs", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeMaxs") })
setMethod("rangeMaxs", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_maxs, FUN.TYPE="numeric")
          })

##' @export rangeWhichMins
setGeneric("rangeWhichMins", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeWhichMins") })
setMethod("rangeWhichMins", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_which_mins, FUN.TYPE="integer")
          })

##' @export rangeWhichMaxs
setGeneric("rangeWhichMaxs", function(x, ir, na.rm=TRUE, simplify=TRUE) { standardGeneric("rangeWhichMaxs") })
setMethod("rangeWhichMaxs", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=TRUE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_which_maxs, FUN.TYPE="integer")
          })


### Internal methods to get directly to summary functions, skipping trim
.rle_view_sums <- function(x, na.rm) { .Call2("RleViews_viewSums", x, na.rm, PACKAGE = "IRanges") }
.rle_view_means <- function(x, na.rm) { .Call2("RleViews_viewMeans", x, na.rm, PACKAGE = "IRanges") }
.rle_view_mins <- function(x, na.rm) { .Call2("RleViews_viewMins", x, na.rm, PACKAGE = "IRanges") }
.rle_view_maxs <- function(x, na.rm) { .Call2("RleViews_viewMaxs", x, na.rm, PACKAGE = "IRanges") }
.rle_view_which_mins <- function(x, na.rm) { .Call2("RleViews_viewWhichMins", x, na.rm, PACKAGE = "IRanges") }
.rle_view_which_maxs <- function(x, na.rm) { .Call2("RleViews_viewWhichMaxs", x, na.rm, PACKAGE = "IRanges") }
