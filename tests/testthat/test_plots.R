library(testthat)
library(genoset)

test_that("We can make genoPlots", {
    data(genoset)
    expect_true( genoPlot( x=genoset.ds,y=genoset.ds[,1,"lrr"] ) )
    expect_true( genoPlot( genoPos(genoset.ds), genoset.ds[,1,"lrr"], locs=rowRanges(genoset.ds) ) )
    expect_true( genoPlot( 1:10, Rle(c(rep(0,5),rep(3,4),rep(1,1))) ) )
})
