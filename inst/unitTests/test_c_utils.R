library(genoset)
library(RUnit)
#library(inline)

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
#  checkEquals(s, mystart)
#  checkEquals(e, myend)
#  
#}

test_RleViews_viewMeans <- function(){

  rle = Rle(1:5, rep(2, 5))
  ir = IRanges(start=c(3, 4), end=c(3, 9), names=c("GENE1", "GENE2"))
  myview = Views(rle, ir)
  s = start(ir)
  w = width(ir)
  one = .Call("RleViews_viewMeans", myview, TRUE, PACKAGE = "IRanges")
  two = .Call("RleViews_viewMeans2", s, w, as.numeric(runValue(rle)), runLength(rle), TRUE, PACKAGE = "genoset")
  check = cbind(as.data.frame(ranges(rle)), runValue(rle))
  checkEquals(one, two, checkNames=FALSE)
  
  rle2 = Rle(c(1, NA, 3, NA, 5), rep(2,5))
  myview2 = Views(rle2, ir)
  one2 = .Call("RleViews_viewMeans", myview2, TRUE, PACKAGE = "IRanges")
  two2 = .Call("RleViews_viewMeans2", s, w, as.numeric(runValue(rle2)), runLength(rle2), TRUE, PACKAGE = "genoset")
  checkEquals(one2,two2, checkNames=FALSE)
  
  one3 = .Call("RleViews_viewMeans", myview2, FALSE, PACKAGE = "IRanges")
  two3 = .Call("RleViews_viewMeans2", s, w, as.numeric(runValue(rle2)), runLength(rle2), FALSE, PACKAGE = "genoset")
  checkEquals(one3,two3, checkNames=FALSE)
}
