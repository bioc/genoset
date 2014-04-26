#### Tests for RleDataFrame and  additional methods for DataFrame

test_RleDataFrame <- function() {
  foo = new("RleDataFrame", listData=list(A=Rle(1:5, rep(2, 5)), B=Rle(6:10,rep(2, 5))), nrows=10L)
  foo2 = new("RleDataFrame", listData=list(A=Rle(1, 2), B=Rle(6,2)), nrows=2L)
  foo3 = RleDataFrame(A=Rle(1:5, rep(2, 5)), B=Rle(6:10,rep(2, 5)))
  checkEquals( foo[1:2, ], foo2, "Row subset")
  checkEquals( foo, foo3, "Create RleDataFrame with new or with RleDataFrame" )
  checkException(new("RleDataFrame", listData=list(1:10, 1:10), nrows=10L), silent=TRUE)
}

test_rowMeans <- function() {
  foo = new("RleDataFrame", listData=list(A=Rle(1:5, rep(2, 5)), B=Rle(6:10,rep(2, 5))), nrows=10L)
  mat = do.call(cbind, lapply(foo, as.numeric))
  checkEquals(rowMeans(mat), rowMeans(foo))
}

test_colMeans <- function() {
  df.ds = DataFrame( a = Rle(c(5,4,3),c(2,2,2)), b = Rle(c(3,6,9),c(1,1,4)) )
  mat.ds = matrix( c(5,5,4,4,3,3,3,6,9,9,9,9), ncol=2, dimnames=list(NULL,c("a","b")))
  checkEquals( colMeans(df.ds), colMeans(mat.ds) )
}
