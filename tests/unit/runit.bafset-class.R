#################################################
# Tests for functions and methods on BAFSet class
#################################################

test_baf2mbaf <- function() {

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
  
  checkEquals( baf2mbaf( baf.ds, hom.cutoff=0.8                                        ), low.cutoff.mbaf.ds, checkNames=FALSE )
  checkEquals( baf2mbaf( baf.ds, hom.cutoff=0.95                                       ), high.cutoff.mbaf.ds, checkNames=FALSE )
  checkEquals( baf2mbaf( baf.ds, calls=ok.calls, call.pairs=all.call.pairs             ), calls.used.mbaf.ds, checkNames=FALSE )
  checkEquals( baf2mbaf( baf.ds, calls=ok.calls, call.pairs=some.call.pairs, hom.cutoff = 0.95  ), some.calls.used.mbaf.ds, checkNames=FALSE )
  checkException( baf2mbaf( baf.ds, calls=rowsmissing.calls, call.pairs=all.call.pairs    ), silent=TRUE )
  checkException( baf2mbaf( baf.ds, calls=extrarows.calls, call.pairs=all.call.pairs      ), silent=TRUE )
  checkException( baf2mbaf(baf.ds, calls=bad.colnames.calls, call.pairs=all.call.pairs ), silent=TRUE )
  
}

