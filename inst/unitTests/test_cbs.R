#########################################
# Test for segmentation-related functions
#########################################
sample.names = LETTERS[11:13]
probe.names = letters[1:10]
locData.rd = RangedData(ranges=IRanges(start=c(1,3,5,7,4,6,2,4,6,8),width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18")
locData.gr = as(locData.rd,"GRanges")
basic.rle.df = DataFrame(
  K = Rle(c(5.3,2.3,1.20),c(4,2,4)),
  L = Rle(c(1.1,1.4,2.2,3.3,0.5),c(1,3,2,2,2)),
  M = Rle(c(3.3,4.3,4.3,6.3,7.3),c(1,3,1,1,4)),
    row.names=rownames(locData.rd))

basic.segs = list(
  K = data.frame( ID = "K", chrom = factor(c("chr1","chr3","chrX"),levels=names(locData.rd)), loc.start = c(1,1,1), loc.end = c(10,10,10),
    num.mark = c(4,2,4), seg.mean = c(5.3,2.3,1.2), stringsAsFactors=FALSE ),
  L = data.frame( ID = "L", chrom = factor(c("chr1","chr1","chr3","chrX","chrX"),levels=names(locData.rd)), loc.start = c(1,3,1,1,5), loc.end = c(2,10,10,4,10),
    num.mark = c(1,3,2,2,2), seg.mean = c(1.1,1.4,2.2,3.3,0.5), stringsAsFactors=FALSE ),
  M = data.frame( ID = "M", chrom = factor(c("chr1","chr1","chr3","chr3","chrX"),levels=names(locData.rd)), loc.start = c(1,3,4,5,1), loc.end = c(2,10,5,6,10),
    num.mark = c(1,3,1,1,4), seg.mean = c(3.3,4.3,4.3,6.3,7.3), stringsAsFactors=FALSE )
  )

basic.segs.after = list(
  K = data.frame( chrom = factor(c("chr1","chr3","chrX"),levels=names(locData.rd)),
    loc.start = c(1,4,2), loc.end = c(7,6,8), num.mark = c(4,2,4), seg.mean = c(5.3,2.3,1.2), stringsAsFactors=FALSE ),
  L = data.frame( chrom = factor(c("chr1","chr1","chr3","chrX","chrX"),levels=names(locData.rd)),
    loc.start = c(1,3,4,2,6), loc.end = c(1,7,6,4,8), num.mark = c(1,3,2,2,2), seg.mean = c(1.1,1.4,2.2,3.3,0.5), stringsAsFactors=FALSE ),
  M = data.frame( chrom = factor(c("chr1","chr1","chr3","chr3","chrX"),levels=names(locData.rd)),
    loc.start = c(1,3,4,6,2), loc.end = c(1,7,4,6,8), num.mark = c(1,3,1,1,4), seg.mean = c(3.3,4.3,4.3,6.3,7.3), stringsAsFactors=FALSE )
  )

stacked.basic.segs.after = data.frame(
  Sample=c(rep("K",3),rep("L",5),rep("M",5)),
  chrom = factor(c("chr1","chr3","chrX","chr1","chr1","chr3","chrX","chrX","chr1","chr1","chr3","chr3","chrX"),levels=names(locData.rd)),
  loc.start = c(1,4,2,1,3,4,2,6,1,3,4,6,2), loc.end = c(7,6,8,1,7,6,4,8,1,7,4,6,8),
  num.mark = c(4,2,4,1,3,2,2,2,1,3,1,1,4),
  seg.mean = c(5.3,2.3,1.2,1.1,1.4,2.2,3.3,0.5,3.3,4.3,4.3,6.3,7.3),
  row.names=c(paste(rep("K",3),1:3,sep="."),paste(rep("L",5),1:5,sep="."),paste(rep("M",5),1:5,sep=".")),
  stringsAsFactors=FALSE)

basic.rds.after = list(
  K = RangedData( ranges=IRanges(start = c(1,4,2), end = c(7,6,8)),space = factor(c("chr1","chr3","chrX"),levels=names(locData.rd)),  score = c(5.3,2.3,1.2), num.mark = c(4,2,4) ),
  L = RangedData( ranges=IRanges(start = c(1,3,4,2,6), end = c(1,7,6,4,8)), space = factor(c("chr1","chr1","chr3","chrX","chrX"),levels=names(locData.rd)),
    score = c(1.1,1.4,2.2,3.3,0.5), num.mark = c(1,3,2,2,2)),
  M = RangedData( ranges=IRanges(start = c(1,3,4,6,2), end = c(1,7,4,6,8)), space = factor(c("chr1","chr1","chr3","chr3","chrX"),levels=names(locData.rd)),
    score = c(3.3,4.3,4.3,6.3,7.3), num.mark = c(1,3,1,1,4))
  )

test_segs2RangedData <- function() {
  checkEquals( segs2RangedData(basic.segs.after$K), basic.rds.after$K )
  checkEquals( segs2RangedData(basic.segs.after$L), basic.rds.after$L )
  checkEquals( segs2RangedData(basic.segs.after$M), basic.rds.after$M )
}

test_segs2Rle <- function() {
  # With RangedData
  checkEquals( segs2Rle( basic.segs[[1]], locData.rd ), basic.rle.df[[1]], checkNames=FALSE )
  checkEquals( segs2Rle( basic.segs[[2]], locData.rd ), basic.rle.df[[2]], checkNames=FALSE )
  checkEquals( segs2Rle( basic.segs[[3]], locData.rd ), basic.rle.df[[3]], checkNames=FALSE )
  na.df = data.frame( chrom = factor(c("chr1","chr1","chr3","chr3","chrX"),levels=names(locData.rd)),
    loc.start = c(2,5,4,6,2), loc.end = c(3,7,4,6,6), num.mark = c(1,2,1,1,3), seg.mean = c(3.3,4.3,4.3,6.3,7.3), stringsAsFactors=FALSE )
  na.rle = Rle( c(NA,3.3,4.3,6.3,7.3,NA), c(1,1,3,1,3,1) )
  checkEquals( segs2Rle( na.df, locData.rd ), na.rle , checkNames=FALSE )

  # With GRanges
  checkEquals( segs2Rle( basic.segs[[1]], locData.gr ), basic.rle.df[[1]], checkNames=FALSE )
  checkEquals( segs2Rle( basic.segs[[2]], locData.gr ), basic.rle.df[[2]], checkNames=FALSE )
  checkEquals( segs2Rle( basic.segs[[3]], locData.gr ), basic.rle.df[[3]], checkNames=FALSE )
  checkEquals( segs2Rle( na.df, locData.gr ), na.rle , checkNames=FALSE )
}

test_segs2RleDataFrame <- function() {
  checkEquals( segs2RleDataFrame( basic.segs, locData.rd ), basic.rle.df, checkNames=FALSE )
  checkEquals( segs2RleDataFrame( basic.segs, locData.gr ), basic.rle.df, checkNames=FALSE )
}

test_segTable <- function() {
  chr.ind = chrIndices(locData.rd)
  start = start(locData.rd)
  end = end(locData.rd)

  # With RangedData
  checkEquals( segTable( basic.rle.df[["K"]], locData.rd), basic.segs.after[["K"]], checkNames=FALSE )
  checkEquals( segTable( basic.rle.df[["L"]], locData.rd), basic.segs.after[["L"]], checkNames=FALSE )
  checkEquals( segTable( basic.rle.df[["M"]], locData.rd), basic.segs.after[["M"]], checkNames=FALSE )
  checkEquals( segTable( basic.rle.df[["M"]], chr.ind=chr.ind, start=start, end=end), basic.segs.after[["M"]], checkNames=FALSE, "segTable on Rle providing chr.ind, start, end" )
  checkEquals( segTable( basic.rle.df, locData.rd ), basic.segs.after, checkNames=FALSE )
  checkEquals( segTable( basic.rle.df, locData.rd, stack=TRUE ), stacked.basic.segs.after, checkNames=FALSE )

  # With GRanges
  checkEquals( segTable( basic.rle.df[["K"]], locData.gr), basic.segs.after[["K"]], checkNames=FALSE )
  checkEquals( segTable( basic.rle.df[["L"]], locData.gr), basic.segs.after[["L"]], checkNames=FALSE )
  checkEquals( segTable( basic.rle.df[["M"]], locData.gr), basic.segs.after[["M"]], checkNames=FALSE )
  checkEquals( segTable( basic.rle.df[["M"]], chr.ind=chr.ind, start=start, end=end), basic.segs.after[["M"]], checkNames=FALSE, "segTable on Rle providing chr.ind, start, end" )
  checkEquals( segTable( basic.rle.df, locData.gr ), basic.segs.after, checkNames=FALSE )
  checkEquals( segTable( basic.rle.df, locData.gr, stack=TRUE ), stacked.basic.segs.after, checkNames=FALSE )
}

test_segPairTable <- function() {
  cn = Rle(c(3,4,5,6),rep(3,4))
  loh = Rle(c(2,4,6,8,10,12),rep(2,6))
  start = c(9:11,4:9,15:17)
  end = start
  locs.rd = RangedData(IRanges(start=start,end=end),space=c(rep("chr1",3),rep("chr2",6),rep("chr3",3)))
  locs.gr = as(locs.rd,"GRanges")
  chr.ind = chrIndices(locs.rd)

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
#  checkIdentical(segs.gr, segPairTable(cn,loh,chr.ind=chr.ind,start=start,end=end))
  
  segs.df = data.frame(
    chrom=factor(c(rep("chr1",2),rep("chr2",4),rep("chr3",2)),levels=c("chr1","chr2","chr3")),
    loc.start = c( 9,11,4,5,7,9,15,16),
    loc.end   = c(10,11,4,6,8,9,15,17),
    num.mark = c(2,1,1,2,2,1,1,2),
    x  = c(3,3,4,4,5,5, 6,  6),
    y = c(2,4,4,6,8,10,10, 12)
    )
  checkEquals(segs.df, segPairTable(cn,loh,chr.ind=chr.ind,start=start,end=end))
  
  cn.df = DataFrame(a=cn,b=cn+1)
  loh.df = DataFrame(a=loh,b=loh+1)
  stacked.segs.df = do.call(rbind,list(a = segPairTable(cn,loh,chr.ind=chr.ind,start=start,end=end), b = segPairTable(cn+1,loh+1,chr.ind=chr.ind,start=start,end=end)))
  stacked.segs.df = cbind(Sample = rep(c("a","b"),each=8),stacked.segs.df)
  checkEquals(stacked.segs.df, segPairTable(cn.df,loh.df,locs=locs.rd,stack=TRUE))
  checkEquals(stacked.segs.df, segPairTable(cn.df,loh.df,locs=locs.gr,stack=TRUE))
}

test.fixSegNAs <- function() {
  x = Rle(c(1,NA,1,5,4,NA,4,2,NA,3), rep(1,10) )
  x.fixed = Rle(c(1,5,4,2,NA,3), c(3,1,3,1,1,1))
  checkIdentical( fixSegNAs(x), x.fixed, "Easy, no NAs on ends" )

  x = Rle(c(1,NA,1,5,4,NA,4,2,NA), rep(1,9) )
  x.fixed = Rle(c(1,5,4,2), c(3,1,3,2))
  checkIdentical( fixSegNAs(x), x.fixed, "NA at end too" )

  x = Rle(c(NA,1,5,4,NA,4,2,NA), rep(1,8) )
  x.fixed = Rle(c(1,5,4,2), c(2,1,3,2))
  checkIdentical( fixSegNAs(x), x.fixed, "NA at beginning too" )

  x = Rle(c(NA,1,5,4,NA,4,2,NA,2), c(1,1,1,2,3,2,1,4,2) )
  x.fixed = Rle(c(1,5,4,2,NA,2), c(2,1,7,1,4,2))
  checkIdentical( fixSegNAs(x), x.fixed, "Longer NA runs" )
}

test_runCBS <- function() {
  sample.names = paste("a",1:2,sep="")
  probe.names =  paste("p",1:30,sep="")
  ds = matrix(c(c(rep(5,20),rep(3,10)),c(rep(2,10),rep(7,10),rep(9,10))),ncol=2,dimnames=list(probe.names,sample.names))
  ds.with.na = matrix(c(c(rep(5,9),NA,rep(5,10),rep(3,10)),c(rep(2,10),rep(7,10),rep(9,10))),ncol=2,dimnames=list(probe.names,sample.names))
  locs.rd = RangedData(ranges=IRanges(start=c(1:20,1:10),width=1,names=probe.names),space=paste("chr",c(rep(1,20),rep(2,10)),sep=""))
  locs.gr = as(locs.rd,"GRanges")
  
  seg.rle.result = DataFrame( a1 = Rle(c(rep(5,20),rep(3,10))), a2 = Rle(c(rep(2,10),rep(7,10),rep(9,10))), row.names=probe.names )
  seg.list.result = list(
    a1 = data.frame( ID=rep("a1",2), chrom=factor(c("chr1","chr2")), loc.start=c(1,1), loc.end=c(20,10), num.mark=c(20,10), seg.mean=c(5,3), stringsAsFactors=FALSE),
    a2 = data.frame( ID=rep("a2",3), chrom=factor(c("chr1","chr1","chr2")), loc.start=c(1,11,1), loc.end=c(10,20,10), num.mark=c(10,10,10), seg.mean=c(2,7,9), stringsAsFactors=FALSE)
    )
  
  # With RangedData
  checkEquals( runCBS(ds,locs.rd, n.cores=8), seg.rle.result, "Return DF of Rle")
  checkEquals( runCBS(ds.with.na,locs.rd, n.cores=8), seg.rle.result, "Return DF of Rle with some NA in starting data")
  checkEquals( runCBS(ds,locs.rd, n.cores=8, return.segs=TRUE), seg.list.result, "Return seg dfs")
  checkEquals( runCBS(seg.rle.result,locs.rd, n.cores=8), seg.rle.result, "Return seg dfs starting from DF of Rle (like mbaf)")
  checkEquals( runCBS(ds,locs.rd, n.cores=8,alpha=0.01), seg.rle.result, "Runs OK with alpha at 0.01 (requires loading of data from DNAcopy)")

  # With RangedData
  checkEquals( runCBS(ds,locs.gr, n.cores=8), seg.rle.result, "Return DF of Rle")
  checkEquals( runCBS(ds.with.na,locs.gr, n.cores=8), seg.rle.result, "Return DF of Rle with some NA in starting data")
  checkEquals( runCBS(ds,locs.gr, n.cores=8, return.segs=TRUE), seg.list.result, "Return seg dfs")
  checkEquals( runCBS(seg.rle.result,locs.gr, n.cores=8), seg.rle.result, "Return seg dfs starting from DF of Rle (like mbaf)")
  checkEquals( runCBS(ds,locs.gr, n.cores=8,alpha=0.01), seg.rle.result, "Runs OK with alpha at 0.01 (requires loading of data from DNAcopy)")
}
