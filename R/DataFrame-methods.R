##### Additional methods for DataFrame class to allow use of DataTable of Rle vectors as assayData element

setClass("RleDataFrame",
         representation(
           rownames = "characterORNULL",
           nrows = "integer"
           ),
         prototype(rownames = NULL,
                   nrows = 0L,
                   listData = structure(list(),  names = character())),
         contains = c("SimpleRleList", "DataFrame")
)

RleDataFrame <- function(..., row.names=NULL, check.names=TRUE) {
  x = DataFrame(..., row.names=row.names, check.names=check.names)
  return( new("RleDataFrame", listData=x@listData, rownames=rownames(x), nrows=nrow(x)) )
}

setMethod("show", "RleDataFrame",
          function(object) {
            message(sprintf("RleDataFrame with %i rows and %i columns\n", nrow(object), ncol(object)))
            show(object@listData)
          })

setGeneric("rowMeans", function(x, na.rm=TRUE, dims=1L) standardGeneric("colMeans") )
setMethod("rowMeans", signature(x="RleDataFrame"),
          function(x) {
            # This will probably have to be .Call, but try in R first
            # set counter for each col to first runLength - 1
            # set vector to first runValues
            # temp = sum vector of first runValues
            # result[1] = temp
            # for i in 2:nrows(x)
            # for cols with zero counter
            #      temp += (next runValue) - (current runValue)
            #      set counter to next runLength
            #      current runValue = next runValue
            # res[i] = temp
            # decrement all counters
            # }
            # divide by ncols at end
            # NAs make this awful of course
            # return(res)
            ## Can use diff() to get running difference of runValues per Rle rather than doing this in the loop
            ##   This will be pretty big, maybe just do it in the loop looking forward to C version
          })

##' Means of columns
##'
##' Calculate means of columns of a DataFrame as if it were a matrix. Allows colmeans
##' in rangeSampleMeans for DataTable just like a real matrix. I'm sure there
##' is much more clever way to do this using aggregate.
##' 
##' @export
##' @param x DataFrame
##' @param na.rm logical
##' @param dims integer
##' @param ... in generic, for extra args in methods.
##' @examples
##'  df.ds = DataFrame( a = Rle(c(5,4,3),c(2,2,2)), b = Rle(c(3,6,9),c(1,1,4)) )
##'  mat.ds = matrix( c(5,5,4,4,3,3,3,6,9,9,9,9), ncol=2, dimnames=list(NULL,c("a","b")))
##'  \dontrun{ identical( colMeans(df.ds), colMeans(mat.ds) ) }
##' @rdname colMeans
setGeneric("colMeans", function(x, na.rm=TRUE, dims=1L) standardGeneric("colMeans") )
##' @rdname colMeans
setMethod("colMeans", signature(x="DataFrame"), function(x,na.rm=TRUE,dims=1L) { return( sapply(x,mean,na.rm=na.rm) ) } )

# Allow eSet constructor to make featureNames from a DataFrame as if it were a matrix
setMethod("annotatedDataFrameFrom",
          signature(object="DataFrame"),
          Biobase:::annotatedDataFrameFromMatrix)
