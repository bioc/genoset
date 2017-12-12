library(genoset)
library(testthat)
#library(inline)

test_leftBound <- function() {
  #.C("leftBound", 1:10, 4, 10, 2)
}

#test_widthToStartEnd <- function() {
#  w = rep(3, 5)
#  s = numeric(5)
#  e = numeric(5)
#  myend = cumsum(w)
#  mystart = myend - 2
#  f <- cfunction(c(width="integer", start="numeric", end="numeric"), "
#                 widthToStartEnd(INTEGER(width), REAL(start), REAL(end), LENGTH(end));
#                 return(R_NilValue);
#", convention=".Call", includes="/gne/home/phaverty/Projects/R/genoset-package/genoset/src/genoset.h")
#  f(w, s, e, length(s))
#  expect_equal(s, mystart)
#  expect_equal(e, myend)
#
#}

test_that("We can calc viewMeans on RLE", {
  rle = Rle(1:5, rep(2, 5))
  ir = IRanges(start=c(3, 4), end=c(3, 9), names=c("GENE1", "GENE2"))
  myview = Views(rle, ir)
  s = start(ir)
  e = end(ir)
  one = .Call("RleViews_viewMeans", myview, TRUE, PACKAGE = "IRanges")
  two = .Call("rangeMeans_rle", s, e, as.numeric(runValue(rle)), runLength(rle), TRUE, PACKAGE = "genoset")
  expect_equivalent(one, two)

  rle2 = Rle(c(1, NA, 3, NA, 5), rep(2,5))
  myview2 = Views(rle2, ir)
  one2 = .Call("RleViews_viewMeans", myview2, TRUE, PACKAGE = "IRanges")
  two2 = .Call("rangeMeans_rle", s, e, as.numeric(runValue(rle2)), runLength(rle2), TRUE, PACKAGE = "genoset")
  expect_equivalent(one2, two2)

  one3 = .Call("RleViews_viewMeans", myview2, FALSE, PACKAGE = "IRanges")
  two3 = .Call("rangeMeans_rle", s, e, as.numeric(runValue(rle2)), runLength(rle2), FALSE, PACKAGE = "genoset")
  expect_equivalent(one3, two3)
})

test_that("We can calc the number of callable bases in each range", {
  rle = Rle(c(3L, 4L, 1L, 5L, 9L), c(4, 2, 3, 4, 5))
  ir = IRanges(start=c(3, 4, 7, 9), end=c(4, 8, 8, 17), names=c("GENE1", "GENE2", "GENE3", "GENE4"))
  s = start(ir)
  e = end(ir)
  one = c(2L, 3L, 0L, 8L)
  two = .Call("numCallable_rle", s, e, runValue(rle), runLength(rle), 3L, PACKAGE = "genoset")
  expect_equal(one, two, checkNames=FALSE)
  expect_equal(one, numCallable(rle, ir, 3L), checkNames=FALSE)
  expect_equal(one, numCallable(rle, cbind(start(ir), end(ir)), 3L), checkNames=FALSE)
  real.mat = cbind(start(ir), end(ir))
  storage.mode(real.mat) = "double"
  expect_equivalent(one, numCallable(rle, real.mat, 3L))
})
