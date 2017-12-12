library(testthat)
library(genoset)

#################################################
# Tests for utility functions
#################################################

test_that("We can swap baf and mbaf", {

  baf.ds = matrix(
    c(1, 0.8, 0.2, 0.3,  0.75, 0.99, 0.5, 0.1,  0.20, 0.1, 0.4, 0.15),
    nrow=4, ncol=3, dimnames = list(c("A","B","C","D"),c("a","b","c")) )

  low.cutoff.mbaf.ds = matrix(
    c(NA, 0.8, 0.8, 0.7,  0.75, NA, 0.5, NA,  0.80, NA, 0.6, NA),
    nrow=4, ncol=3, dimnames = list(c("A","B","C","D"),c("a","b","c")) )

  high.cutoff.mbaf.ds = matrix(
    c(NA, 0.8, 0.8, 0.7,  0.75, NA, 0.5, 0.9,  0.80, 0.9, 0.6, 0.85),
    nrow=4, ncol=3, dimnames = list(c("A","B","C","D"),c("a","b","c")) )

  calls.used.mbaf.ds = matrix(
    c(NA, 0.8, NA, NA,  NA, NA, 0.5, NA,  0.8, 0.9, NA, NA),
    nrow=4, ncol=3, dimnames = list(c("A","B","C","D"),c("a","b","c")) )

  rowsmissing.calls.used.mbaf.ds = matrix(
    c(NA, 0.8, NA, NA,  NA, NA, 0.5, NA,  0.8, 0.9, NA, NA),
    nrow=4, ncol=3, dimnames = list(c("A","B","C","D"),c("a","b","c")) )

  some.calls.used.mbaf.ds = matrix(
    c(NA, 0.8, NA, NA,  0.75, NA, 0.5, 0.9,  0.8, 0.9, NA, NA),
    nrow=4, ncol=3, dimnames = list(c("A","B","C","D"),c("a","b","c")) )

  all.call.pairs = c("a","b","c")
  names(all.call.pairs) = c("a","b","c")

  some.call.pairs = c("a","c")
  names(some.call.pairs) = c("a","c")

  ok.calls = matrix(
    c("AA","AT","GG","TT",  "AA","GC","TA","AA", "AT","AT","AA","TT"),
    nrow=4, ncol=3, dimnames = list(c("A","B","C","D"),c("a","b","c")) )

  rowsmissing.calls = matrix(
    c("CC","CC","GG",  "TT","AT","TT", "TA","CC","AA"),
    nrow=3, ncol=3, dimnames = list(c("A","C","D"),c("a","b","c")) )

  extrarows.calls = matrix(
    c("TT","AT","CC","AA","CC",  "AA","TA","AT","CC","AA", "GC","AT","CC","AA","TA"),
    nrow=5, ncol=3, dimnames = list(c("A","B","C","D","E"),c("a","b","c")) )

  bad.colnames.calls = matrix(
    c("TT","AT","CC","AG",  "CC","AG","CG","AA", "AT","AC","CC","TT"),
    nrow=4, ncol=3, dimnames=list(c("A","B","C","D"),c("a","b","f")) )

  expect_equal( baf2mbaf( baf.ds, hom.cutoff=0.8                                        ), low.cutoff.mbaf.ds, checkNames=FALSE )
  expect_equal( baf2mbaf( baf.ds, hom.cutoff=0.95                                       ), high.cutoff.mbaf.ds, checkNames=FALSE )
  expect_equal( baf2mbaf( baf.ds, calls=ok.calls, call.pairs=all.call.pairs             ), calls.used.mbaf.ds, checkNames=FALSE )
  expect_equal( baf2mbaf( baf.ds, calls=ok.calls, call.pairs=some.call.pairs, hom.cutoff = 0.95  ), some.calls.used.mbaf.ds, checkNames=FALSE )
  expect_error( baf2mbaf( baf.ds, calls=rowsmissing.calls, call.pairs=all.call.pairs    ), silent=TRUE )
  expect_error( baf2mbaf( baf.ds, calls=extrarows.calls, call.pairs=all.call.pairs      ), silent=TRUE )
  expect_error( baf2mbaf(baf.ds, calls=bad.colnames.calls, call.pairs=all.call.pairs ), silent=TRUE )

})

test_that("We can correct for GC", {

  input.vector = c(rep(0.05,50),rep(0.08,50))
  gc = input.vector
  output.vector = rep(0,100)
  expect_equal( gcCorrect(input.vector, gc, retain.mean=FALSE ), output.vector )
  expect_equal( gcCorrect(input.vector, gc, retain.mean=TRUE ), output.vector + mean(input.vector) )

  input.matrix = matrix(c(input.vector,input.vector),ncol=2)
  output.matrix = matrix(c(output.vector,output.vector),ncol=2)
  expect_equal( gcCorrect(input.matrix, gc, retain.mean=FALSE ), output.matrix )

  input.matrix.w.na = input.matrix
  input.matrix.w.na[ c(25,75),  ] = NA
  output.matrix.w.na = output.matrix
  output.matrix.w.na[ c(25,75), ] = NA
  expect_equal( gcCorrect(input.matrix.w.na, gc, retain.mean=FALSE ), output.matrix.w.na )
})

test_that("We can swap cn and lr", {
  ploidy = 2:3
  dimnames = list(letters[1:2], LETTERS[1:2])
  cn = matrix(2:5, byrow=TRUE, ncol=2, dimnames=dimnames)
  cn[1, 2] = -0.1
  lr = log2( matrix(2:5, byrow=TRUE, ncol=2, dimnames=dimnames) / 2 )
  lr[1, 2] = log2(0.001 / 2)
  relative.lr = matrix(c(0, log2(0.001/3), log2(4/2), log2(5/3)), byrow=TRUE, ncol=2, dimnames=dimnames)
  expect_equal( lr, cn2lr(cn), label = "Assume diploid")
  expect_equal( relative.lr, cn2lr(cn, ploidy), label = "Use ploidy")
  expect_equal( relative.lr[, 2], cn2lr(cn[, 2], ploidy[2]), label = "Use ploidy one sample")
  expect_error( cn2lr(cn, 1:8), label = "Ploidy and cn must match in size", silent=TRUE)

  cn.df = as(cn, "DataFrame")
  lr.df = as(lr, "DataFrame")
  relative.lr.df = as(relative.lr, "DataFrame")

#  expect_equal( lr.df, cn2lr(cn.df), "Assume diploid")
#  expect_equal( relative.lr.df, cn2lr(cn.df, ploidy), "Use ploidy")
  expect_error( cn2lr(cn.df, 1:8), label = "Ploidy and cn must match in size", silent=TRUE)
  expect_equal( lr2cn(c(-1, 0, 1)), c(1, 2, 4) )
})

test_that("We can mode center", {
    x = matrix( c(1, 1, 2, 2, 2, 2, 3, 2, 1), ncol=3 )
    x2 = matrix( c(0, 0, 1, 0, 0, 0, 1, 0, -1), ncol=3 )
    expect_equal( modeCenter( x ), x2, tolerance=0.05 )
})

test_that("We can rbind data frames", {
    df = data.frame(a=letters,b=LETTERS,c=1:length(letters),d=factor(letters), stringAsFactors=FALSE)
    df_list = lapply(paste0(letters,LETTERS), function(x) { return(df) })
    expect_equal( rbindDataframe(df_list), do.call(rbind,df_list) )
})
