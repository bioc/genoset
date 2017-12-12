#########################################
# Test for segmentation-related functions
#########################################
library(testthat)
library(genoset)

sample.names = LETTERS[11:13]
probe.names = letters[1:10]
locData.gr = GRanges(ranges=IRanges(start=c(1,3,5,7,4,6,2,4,6,8),width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
basic.rle.df = RleDataFrame(
  K = Rle(c(5.3,2.3,1.20),c(4,2,4)),
  L = Rle(c(1.1,1.4,2.2,3.3,0.5),c(1,3,2,2,2)),
  M = Rle(c(3.3,4.3,4.3,6.3,7.3),c(1,3,1,1,4)),
    row.names=names(locData.gr))

basic.segs = list(
  K = data.frame( ID = "K", chrom = factor(c("chr1","chr3","chrX"),levels=chrNames(locData.gr)), loc.start = c(1,1,1), loc.end = c(10,10,10),
    num.mark = c(4,2,4), seg.mean = c(5.3,2.3,1.2), stringsAsFactors=FALSE ),
  L = data.frame( ID = "L", chrom = factor(c("chr1","chr1","chr3","chrX","chrX"),levels=chrNames(locData.gr)), loc.start = c(1,3,1,1,5), loc.end = c(2,10,10,4,10),
    num.mark = c(1,3,2,2,2), seg.mean = c(1.1,1.4,2.2,3.3,0.5), stringsAsFactors=FALSE ),
  M = data.frame( ID = "M", chrom = factor(c("chr1","chr1","chr3","chr3","chrX"),levels=chrNames(locData.gr)), loc.start = c(1,3,4,5,1), loc.end = c(2,10,5,6,10),
    num.mark = c(1,3,1,1,4), seg.mean = c(3.3,4.3,4.3,6.3,7.3), stringsAsFactors=FALSE )
  )

basic.segs.after = list(
  K = data.frame( chrom = factor(c("chr1","chr3","chrX"),levels=chrNames(locData.gr)),
    loc.start = c(1,4,2), loc.end = c(7,6,8), num.mark = c(4,2,4), seg.mean = c(5.3,2.3,1.2), stringsAsFactors=FALSE ),
  L = data.frame( chrom = factor(c("chr1","chr1","chr3","chrX","chrX"),levels=chrNames(locData.gr)),
    loc.start = c(1,3,4,2,6), loc.end = c(1,7,6,4,8), num.mark = c(1,3,2,2,2), seg.mean = c(1.1,1.4,2.2,3.3,0.5), stringsAsFactors=FALSE ),
  M = data.frame( chrom = factor(c("chr1","chr1","chr3","chr3","chrX"),levels=chrNames(locData.gr)),
    loc.start = c(1,3,4,6,2), loc.end = c(1,7,4,6,8), num.mark = c(1,3,1,1,4), seg.mean = c(3.3,4.3,4.3,6.3,7.3), stringsAsFactors=FALSE )
  )

stacked.basic.segs.after = data.frame(
  chrom = factor(c("chr1","chr3","chrX","chr1","chr1","chr3","chrX","chrX","chr1","chr1","chr3","chr3","chrX"),levels=chrNames(locData.gr)),
  loc.start = c(1,4,2,1,3,4,2,6,1,3,4,6,2), loc.end = c(7,6,8,1,7,6,4,8,1,7,4,6,8),
  num.mark = c(4,2,4,1,3,2,2,2,1,3,1,1,4),
  seg.mean = c(5.3,2.3,1.2,1.1,1.4,2.2,3.3,0.5,3.3,4.3,4.3,6.3,7.3),
  Sample=factor(c(rep("K",3),rep("L",5),rep("M",5))),
  row.names=NULL,
  stringsAsFactors=FALSE)

basic.rds.after = list(
  K = GRanges( ranges=IRanges(start = c(1,4,2), end = c(7,6,8)),seqnames = factor(c("chr1","chr3","chrX"),levels=chrNames(locData.gr)), num.mark = c(4,2,4), seg.mean = c(5.3,2.3,1.2) ),
  L = GRanges( ranges=IRanges(start = c(1,3,4,2,6), end = c(1,7,6,4,8)), seqnames = factor(c("chr1","chr1","chr3","chrX","chrX"),levels=chrNames(locData.gr)),
    num.mark = c(1,3,2,2,2), seg.mean = c(1.1,1.4,2.2,3.3,0.5)),
  M = GRanges( ranges=IRanges(start = c(1,3,4,6,2), end = c(1,7,4,6,8)), seqnames = factor(c("chr1","chr1","chr3","chr3","chrX"),levels=chrNames(locData.gr)),
    num.mark = c(1,3,1,1,4), seg.mean = c(3.3,4.3,4.3,6.3,7.3))
  )

test_that("We can convert segs to GRanges", {
  expect_equal( segs2Granges(basic.segs.after$K), basic.rds.after$K )
  expect_equal( segs2Granges(basic.segs.after$L), basic.rds.after$L )
  expect_equal( segs2Granges(basic.segs.after$M), basic.rds.after$M )
})

test_that("We can convert segs to Rle", {
  # With GRanges
  na.df = data.frame( chrom = factor(c("chr1", "chr1", "chr3", "chr3", "chrX"), levels=chrNames(locData.gr)),
    loc.start = c(3, 5, 4, 6, 2),  loc.end = c(3, 7, 4, 6, 6),  num.mark = c(1, 2, 1, 1, 3),  seg.mean = c(3.3, 4.3, 4.3, 6.3, 7.3),  stringsAsFactors=FALSE )
  na.rle = Rle( c(NA, 3.3, 4.3, 6.3, 7.3, NA),  c(1, 1, 3, 1, 3, 1) )
  expect_equal( segs2Rle( basic.segs[[1]], locData.gr ), basic.rle.df[[1]], checkNames=FALSE )
  expect_equal( segs2Rle( basic.segs[[2]], locData.gr ), basic.rle.df[[2]], checkNames=FALSE )
  expect_equal( segs2Rle( basic.segs[[3]], locData.gr ), basic.rle.df[[3]], checkNames=FALSE )
  expect_equal( segs2Rle( na.df, locData.gr ), na.rle , checkNames=FALSE )
})

test_that("We can convert segs to DF", {
  expect_equal( segs2RleDataFrame( basic.segs, locData.gr ), basic.rle.df, checkNames=FALSE )
})

test_segTable <- function() {
  chr.ind = chrIndices(locData.gr)
  start = start(locData.gr)
  end = end(locData.gr)

  # With GRanges
  expect_equal( segTable( basic.rle.df[["K"]], locData.gr), basic.segs.after[["K"]], checkNames=FALSE )
  expect_equal( segTable( basic.rle.df[["L"]], locData.gr), basic.segs.after[["L"]], checkNames=FALSE )
  expect_equal( segTable( basic.rle.df[["M"]], locData.gr), basic.segs.after[["M"]], checkNames=FALSE )
  expect_equal( segTable( basic.rle.df[["M"]], chr.ind=chr.ind, start=start, end=end), basic.segs.after[["M"]], checkNames=FALSE, label = "segTable on Rle providing chr.ind, start, end" )
  expect_equal( segTable( basic.rle.df, locData.gr ), basic.segs.after, checkNames=FALSE )
  expect_equal( segTable( basic.rle.df, locData.gr, stack=TRUE ), stacked.basic.segs.after, checkNames=FALSE )
}

test_that("We can convert segs to pair table", {
  cn = Rle(c(3,4,5,6),rep(3,4))
  loh = Rle(c(2,4,6,8,10,12),rep(2,6))
  start = c(9:11,4:9,15:17)
  end = start
  locs.gr = GRanges(IRanges(start=start,end=end),seqnames=c(rep("chr1",3),rep("chr2",6),rep("chr3",3)))
  chr.ind = chrIndices(locs.gr)

  # Test returning a GRanges from segPairTable
#  segs.gr = GRanges(
#    ranges = IRanges(
#      start = c(9, 11,4,5,7,9,15,16),
#      end   = c(10,11,4,6,7,9,15,17)
#      ),
#    seqnames=c(rep("chr1",2),rep("chr2",4),rep("chr3",2)),
##    seqlengths=list("chr1"=11,"chr2"=9,"chr3"=17),
#    x  = c(3,3,4,4,5,5, 6,  6),
#    y = c(2,4,4,6,8,10,10, 12)
#    )
#  expect_identical(segs.gr, segPairTable(cn,loh,chr.ind=chr.ind,start=start,end=end))

  segs.df = data.frame(
    chrom=factor(c(rep("chr1",2),rep("chr2",4),rep("chr3",2)),levels=c("chr1","chr2","chr3")),
    loc.start = c( 9,11,4,5,7,9,15,16),
    loc.end   = c(10,11,4,6,8,9,15,17),
    num.mark = c(2,1,1,2,2,1,1,2),
    x  = c(3,3,4,4,5,5, 6,  6),
    y = c(2,4,4,6,8,10,10, 12)
    )
  expect_equal(segs.df, segPairTable(cn,loh,chr.ind=chr.ind,start=start,end=end))

  cn.df = DataFrame(a=cn,b=cn+1)
  loh.df = DataFrame(a=loh,b=loh+1)
  stacked.segs.df = do.call(rbind,list(a = segPairTable(cn,loh,chr.ind=chr.ind,start=start,end=end), b = segPairTable(cn+1,loh+1,chr.ind=chr.ind,start=start,end=end)))
  stacked.segs.df = cbind(
    stacked.segs.df,
    Sample = factor(rep(c("a","b"),each=8)),
    stringsAsFactors=FALSE,row.names=NULL)
  expect_equal(stacked.segs.df, segPairTable(cn.df,loh.df,locs=locs.gr,stack=TRUE))
})

test_that("We can fix NA segments", {
  x = Rle(c(1,NA,1,5,4,NA,4,2,NA,3), rep(1,10) )
  x.fixed = Rle(c(1,5,4,2,NA,3), c(3,1,3,1,1,1))
  expect_identical( fixSegNAs(x), x.fixed, "Easy, no NAs on ends" )

  x = Rle(c(1,NA,1,5,4,NA,4,2,NA), rep(1,9) )
  x.fixed = Rle(c(1,5,4,2), c(3,1,3,2))
  expect_identical( fixSegNAs(x), x.fixed, "NA at end too" )

  x = Rle(c(NA,1,5,4,NA,4,2,NA), rep(1,8) )
  x.fixed = Rle(c(1,5,4,2), c(2,1,3,2))
  expect_identical( fixSegNAs(x), x.fixed, "NA at beginning too" )

  x = Rle(c(NA,1,5,4,NA,4,2,NA,2), c(1,1,1,2,3,2,1,4,2) )
  x.fixed = Rle(c(1,5,4,2,NA,2), c(2,1,7,1,4,2))
  expect_identical( fixSegNAs(x), x.fixed, "Longer NA runs" )
})

test_that("We can run CBS", {
  sample.names = paste("a",1:2,sep="")
  probe.names =  paste("p",1:30,sep="")
  ds = matrix(c(c(rep(5,20),rep(3,10)),c(rep(2,10),rep(7,10),rep(9,10))),ncol=2,dimnames=list(probe.names,sample.names))
  ds.with.na = matrix(c(c(rep(5,9),NA,rep(5,10),rep(3,10)),c(rep(2,10),rep(7,10),rep(9,10))),ncol=2,dimnames=list(probe.names,sample.names))
  locs.gr = GRanges(ranges=IRanges(start=c(1:20,1:10),width=1,names=probe.names),seqnames=paste("chr",c(rep(1,20),rep(2,10)),sep=""))

  seg.rle.result = RleDataFrame( a1 = Rle(c(rep(5,20),rep(3,10))), a2 = Rle(c(rep(2,10),rep(7,10),rep(9,10))), row.names=probe.names )
  seg.list.result = list(
    a1 = data.frame( ID=rep("a1",2), chrom=factor(c("chr1","chr2")), loc.start=c(1,1), loc.end=c(20,10), num.mark=c(20,10), seg.mean=c(5,3), stringsAsFactors=FALSE),
    a2 = data.frame( ID=rep("a2",3), chrom=factor(c("chr1","chr1","chr2")), loc.start=c(1,11,1), loc.end=c(10,20,10), num.mark=c(10,10,10), seg.mean=c(2,7,9), stringsAsFactors=FALSE)
    )

  # With GRanges
  expect_equal( runCBS(ds,locs.gr, n.cores=8), seg.rle.result, label = "Return DF of Rle")
  expect_equal( runCBS(ds.with.na,locs.gr, n.cores=8), seg.rle.result, label = "Return DF of Rle with some NA in starting data")
  expect_equal( runCBS(ds,locs.gr, n.cores=8, return.segs=TRUE), seg.list.result, label = "Return seg dfs")
  expect_equal( runCBS(seg.rle.result,locs.gr, n.cores=8), seg.rle.result, label = "Return seg dfs starting from DF of Rle (like mbaf)")
  expect_equal( runCBS(ds,locs.gr, n.cores=8,alpha=0.01), seg.rle.result, label = "Runs OK with alpha at 0.01 (requires loading of data from DNAcopy)")
})

test_that("We can convert segs to GRanges", {
  segs = data.frame(loc.start=1:4, loc.end=c(5, 7, 9, 22), chrom=factor(c("1", "1", "2", "2"), levels=c("1", "2")), num.mark=letters[1:4], goo=LETTERS[1:4], stringsAsFactors=FALSE)
  gr = GRanges(IRanges(start=1:4, end=c(5, 7, 9, 22)), seqnames=factor(c("1", "1", "2", "2"), levels=c("1", "2")), num.mark=letters[1:4], goo=LETTERS[1:4])
  expect_identical(gr, segs2Granges(segs))
})

test_that("We can figure out average segment lengths", {
  segs = data.frame(loc.start=1:4, loc.end=c(5, 7, 9, 22), chrom=factor(c("1", "1", "2", "2"), levels=c("1", "2")), num.mark=letters[1:4], goo=LETTERS[1:4], stringsAsFactors=FALSE)
  seg.list = list(
    foo=data.frame(loc.start=1:4, loc.end=c(5, 7, 9, 22), chrom=factor(c("1", "1", "2", "2"), levels=c("1", "2")), num.mark=letters[1:4], goo=LETTERS[1:4], stringsAsFactors=FALSE),
    goo=data.frame(loc.start=2:5, loc.end=c(5, 7, 9, 20), chrom=factor(c("1", "1", "2", "2"), levels=c("1", "2")), num.mark=letters[1:4], goo=LETTERS[1:4], stringsAsFactors=FALSE)
    )
  gene.gr = GRanges(IRanges(start=c(1, 3), end=c(1, 20), names=c("EGFR", "ERBB2")), seqnames=factor(c("1", "2"), levels=c("1", "2")))
  len = c("EGFR"=5, "ERBB2"=13)
  len.mat = matrix(c(5, 13, 4, 11), ncol=2, dimnames=list(c("EGFR", "ERBB2"), c("foo", "goo")))
  expect_identical(len, rangeSegMeanLength(gene.gr, segs))
  expect_identical(len.mat, rangeSegMeanLength(gene.gr, seg.list))
})
