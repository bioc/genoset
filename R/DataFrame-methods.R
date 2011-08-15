##### Additional methods for DataFrame class to allow use of DataTable of Rle vectors as assayData element


# Allow colmeans in rangeSampleMeans for DataTable just like a real matrix
# I'm sure there is much more clever way to do this using window or something

##' Means of columns
##'
##' Get means of columns of a DataFrame as if it were a matrix
##' 
##' @export colMeans
##' @param x DataFrame
##' @param na.rm logical
##' @param dims integer
##' @author Peter M. Haverty
##' @examples
##'  df.ds = DataFrame( a = Rle(c(5,4,3),c(2,2,2)), b = Rle(c(3,6,9),c(1,1,4)) )
##'  mat.ds = matrix( c(5,5,4,4,3,3,3,6,9,9,9,9), ncol=2, dimnames=list(NULL,c("a","b")))
##'  identical( colMeans(df.ds), colMeans(mat.ds) )
##' @rdname colMeans
setGeneric("colMeans", function(x,na.rm=TRUE,dims=1L) standardGeneric("colMeans") )
##' @rdname colMeans
setMethod("colMeans", signature(x="DataFrame"), function(x,na.rm=TRUE,dims=1L) { return( sapply(x,mean,na.rm=na.rm) ) } )

# Allow eSet constructor to make featureNames from a DataFrame as if it were a matrix
setMethod("annotatedDataFrameFrom",
          signature(object="DataFrame"),
          Biobase:::annotatedDataFrameFromMatrix)

# Allow eSet constructor to make featureNames from a big.matrix as if it were a matrix
setMethod("annotatedDataFrameFrom",
          signature(object="big.matrix"),
          Biobase:::annotatedDataFrameFromMatrix)

