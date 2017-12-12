#### Tests for RleDataFrame and  additional methods for DataFrame
library(genoset)
library(testthat)

test_that("We can make RLEDFs", {
  foo = new("RleDataFrame", listData=list(A=Rle(1:5, rep(2, 5)), B=Rle(6:10,rep(2, 5))), nrows=10L)
  foo2 = new("RleDataFrame", listData=list(A=Rle(1, 2), B=Rle(6,2)), nrows=2L)
  foo3 = RleDataFrame(A=Rle(1:5, rep(2, 5)), B=Rle(6:10,rep(2, 5)))
  foo4 = DataFrame(A=Rle(1:5, rep(2, 5)), B=Rle(6:10,rep(2, 5)))
  expect_equal( foo3, as(foo3, "RleDataFrame"), label = "Coercion from DataFrame to RleDataFrame")
  expect_equal( RleList(foo3@listData, compress=FALSE), as(foo3, "RleList"), label = "Coercion to a List")
  expect_equal( as.list(foo3), foo3@listData, label = "Coercion from RleDataFrame to list")
  expect_equal( as.data.frame(foo3), data.frame(lapply(foo3, as.vector)), label = "data.frame coercion")
  expect_equal( as.matrix(as.data.frame(foo3)), as(foo3, "matrix"), label = "Matrix coercion")
  expect_equal( as.matrix(as.data.frame(foo3)), as.matrix(foo3), label = "Matrix coercion again")
  expect_equal( data.frame(A=foo3[[1]], B=foo3[[2]]), as.data.frame(foo3), label = "Coercion to data.frame")
  expect_equal( foo[1:2, ], foo2, label = "Row subset")
  expect_equal( foo, foo3, label = "Create RleDataFrame with new or with RleDataFrame" )
  expect_error(new("RleDataFrame", listData=list(1:10, 1:10), nrows=10L), silent=TRUE)
})

test_that("We can summarize rows", {
  foo = new("RleDataFrame", listData=list(A=Rle(c(NA, 2:3, NA, 5), rep(2, 5)), B=Rle(c(6:7, NA, 8:10),c(3,2,1,2,1,1))), nrows=10L)
  mat = do.call(cbind, lapply(foo, as.numeric))
  expect_equal(rowMeans(mat), as.numeric(rowMeans(foo)))
  expect_equal(rowMeans(mat, na.rm=TRUE), as.numeric(rowMeans(foo, na.rm=TRUE)))
  expect_equal(rowSums(mat), as.numeric(rowSums(foo)))
  expect_equal(rowSums(mat, na.rm=TRUE), as.numeric(rowSums(foo, na.rm=TRUE)))
})

test_that("We can summarize cols", {
  df.ds = DataFrame( a = Rle(c(5,4,3),c(2,2,2)), b = Rle(c(3,6,9),c(1,1,4)) )
  rle.df = new("RleDataFrame", listData=df.ds@listData, nrows=nrow(df.ds))
  mat.ds = matrix( c(5,5,4,4,3,3,3,6,9,9,9,9), ncol=2, dimnames=list(NULL,c("a","b")))
  expect_equal( suppressWarnings(colMeans(df.ds)), colMeans(mat.ds) )
  expect_equal( colMeans(rle.df), colMeans(mat.ds) )
})
