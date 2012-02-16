#### Tests for additional methods for DataFrame

test_colMeans <- function() {
  df.ds = DataFrame( a = Rle(c(5,4,3),c(2,2,2)), b = Rle(c(3,6,9),c(1,1,4)) )
  mat.ds = matrix( c(5,5,4,4,3,3,3,6,9,9,9,9), ncol=2, dimnames=list(NULL,c("a","b")))
  checkEquals( colMeans(df.ds), colMeans(mat.ds) )
}
