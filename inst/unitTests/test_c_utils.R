#library(inline)
#
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
    one = .Call2("RleViews_viewMeans", myview, TRUE, PACKAGE = "IRanges")
#    two = .Call2("RleViews_viewMeans", s, w, as.numeric(runValue(rle)), runLength(rle), TRUE, PACKAGE = "genoset")
    three = .Call2("RleViews_viewMeans2", s, w, as.numeric(runValue(rle)), runLength(rle), TRUE, PACKAGE = "genoset")
    check = cbind(as.data.frame(ranges(rle)), runValue(rle))
#    checkEquals(one, two)
    checkEquals(one, three, checkNames=FALSE)
    return(TRUE)
  }
