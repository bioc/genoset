### Methods to do view operations on an RleDataFrame
###  viewMeans(Views(rledf,iranges)) works, but gives list.  Consistency good, but I really wanna simplify via vapply here ...

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
           }, USE.NAMES=TRUE, 
           FUN.VALUE=structure( vector( FUN.TYPE, nviews ), names=names(ir) ) )
  } else {
    val = lapply(x,   	 
      function(rle) {
        RLEFUN(myviewfun(rle, ir), na.rm=na.rm)
      })
    }
  return(val)
}

##' Calculate min/max/sum/mean/whichmin/whichmax over ranges on each column of an RleDataFrame.
##'
##' Loop over the Rle objects in an RleDataFrame, calculate the appropriate statistic for each view. If simplify == FALSE,
##' this function returns a vector for each Rle. If simplify == TRUE, it returns a vector for the case of a single range, otherwise,
##' a matrix. Rownames for the matrix are taken from the names of the argument \code{ir}.
##' @param x RleDataFrame
##' @param ir IRanges or matrix, views on every Rle. If \code{ir} is a matrix, the first
##' two columns are used as as the starts and stops. Names for the ranges are taken from rownames of the matrix. Such a matrix can be
##' constructed with \code{boundingIndicesByChr}.
##' @param na.rm scalar logical, ignore NAs in calculations?
##' @param simplify scalar logical, simplify result? For a single view, a vector, otherwise a matrix with one row per view.
##' @param RLEFUN function, internal rle view summary function like .rle_range_means
##' @param FUN.TYPE scalar character, the storage mode for the returned vector or matrix (when simplify==TRUE).
##' @return With \code{simplify == TRUE}, a vector for single view or a matrix
##' otherwise. When \code{simplify == FALSE}, a list of vectors length ncol(x) where each element is of length \code{nrows(ir)}.
##' @keywords internal
##' @rdname do_rledf_range_summary
##' @seealso RleDataFrame boundingIndicesByChr
##' @family views
.do_rledf_range_summary <- function(x, ir, na.rm=FALSE, simplify=TRUE, RLEFUN, FUN.TYPE=c("numeric", "double", "integer", "logical")) {
  # Make an IRanges from ranges matrix if necessary
  if (is.matrix(ir)) {
    start=ir[, 1]
    end=ir[, 2]
    names=rownames(ir)
  } else if (is(ir, "IRanges")) {
    start = start(ir)
    end = end(ir)
    names = names(ir)
  } else {
    stop("x must be a two-column matrix or an IRanges.")
  }
  # Trim IRanges once if necessary
  start[ start < 1L ] = 1L
  end[ end > nrow(x) ] = nrow(x)
  width = (end - start) + 1L
  # Calculate the view stats
  if (simplify == TRUE) {
    FUN.TYPE = match.arg(FUN.TYPE)
    val = vapply(x,   	 
           FUN=function(rle) {
             RLEFUN(start, width, runValue(rle), runLength(rle), na.rm=na.rm)
           }, USE.NAMES=TRUE, 
           FUN.VALUE=structure( vector( FUN.TYPE, length(start) ), names=names) )
  } else {
    val = lapply(x,   	 
      function(rle) {
        structure( RLEFUN(start, width, runValue(rle), runLength(rle), na.rm=na.rm), names=names)
      })
    }
  return(val)
}

##' @export rangeSums
setGeneric("rangeSums", function(x, ir, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeSums") })
setMethod("rangeSums", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_sums, FUN.TYPE="numeric")
          })

##' @export rangeMins
setGeneric("rangeMins", function(x, ir, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeMins") })
setMethod("rangeMins", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_mins, FUN.TYPE="numeric")
          })

##' @export rangeMaxs
setGeneric("rangeMaxs", function(x, ir, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeMaxs") })
setMethod("rangeMaxs", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_maxs, FUN.TYPE="numeric")
          })

##' @export rangeWhichMins
setGeneric("rangeWhichMins", function(x, ir, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeWhichMins") })
setMethod("rangeWhichMins", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_which_mins, FUN.TYPE="integer")
          })

##' @export rangeWhichMaxs
setGeneric("rangeWhichMaxs", function(x, ir, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeWhichMaxs") })
setMethod("rangeWhichMaxs", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_views(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_view_which_maxs, FUN.TYPE="integer")
          })


##' @export rangeMeans
setGeneric("rangeMeans", function(x, ir, na.rm=FALSE, simplify=TRUE) { standardGeneric("rangeMeans") })
setMethod("rangeMeans", signature=signature(x="RleDataFrame"), 
          function(x, ir, na.rm=FALSE, simplify=TRUE) {
            .do_rledf_range_summary(x, ir, na.rm=na.rm, simplify=simplify, RLEFUN=.rle_range_means, FUN.TYPE="numeric")
          })

setMethod("rangeMeans", signature=signature(x="vector"), 
          function(x, bounds, na.rm=FALSE) {
              if (! is.matrix(bounds) && ncol(bounds) == 2) {
                  stop("bounds must be a matrix with 2 columns\n")
              }
              if (!is.double(x)) {
                  storage.mode(x) = "double"
              }
              if (!is.integer(bounds)) {
                  storage.mode(bounds) = "integer"
              }
              ans = .Call("rangeMeans_vector", x, bounds)
              return(ans)
          })

setMethod("rangeMeans", signature=signature(x="ANY"),
          function(x, all.indices, na.rm=FALSE) {
              range.means = vapply( structure(seq.int(length.out=ncol(data.matrix)), names=colnames(data.matrix)),
                  FUN=function(x) { rangeMeans(as.numeric(data.matrix[, x])) },
                  FUN.VALUE = structure(numeric(nrow(data.matrix)), names=rownames(all.indices)) )
          return(range.means)
      })

##' @export rangeColMeans
rangeColMeans <- function(x, all.indices) {
.Deprecated("rangeMeans", "rangeColMeans has changed to rangeMeans.")
rangeMeans(x, all.indices)
}

### Internal methods to get directly to summary functions, using Views, but skipping trim
.rle_view_sums <- function(x, na.rm) { .Call("RleViews_viewSums", x, na.rm, PACKAGE = "IRanges") }
.rle_view_means <- function(x, na.rm) { .Call("RleViews_viewMeans", x, na.rm, PACKAGE = "IRanges") }
.rle_view_mins <- function(x, na.rm) { .Call("RleViews_viewMins", x, na.rm, PACKAGE = "IRanges") }
.rle_view_maxs <- function(x, na.rm) { .Call("RleViews_viewMaxs", x, na.rm, PACKAGE = "IRanges") }
.rle_view_which_mins <- function(x, na.rm) { .Call("RleViews_viewWhichMins", x, na.rm, PACKAGE = "IRanges") }
.rle_view_which_maxs <- function(x, na.rm) { .Call("RleViews_viewWhichMaxs", x, na.rm, PACKAGE = "IRanges") }

### Internal methods to get directly to summary functions, skipping trim and Views
.rle_range_means <- function(start, width, values, lengths, na.rm) {.Call("RleViews_viewMeans2", start, width, as.numeric(values), lengths, na.rm=na.rm, PACKAGE = "genoset") }
