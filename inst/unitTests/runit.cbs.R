#########################################
# Test for segmentation-related functions
#########################################
sample.names = LETTERS[11:13]
probe.names = letters[1:10]
locData.rd = RangedData(ranges=IRanges(start=c(1,3,5,7,4,6,2,4,6,8),width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18")
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
  chrom = factor(c("chr1","chr3","chrX","chr1","chr1","chr3","chrX","chrX","chr1","chr1","chr3","chr3","chrX"),levels=names(locData.rd)),
  loc.start = c(1,4,2,1,3,4,2,6,1,3,4,6,2), loc.end = c(7,6,8,1,7,6,4,8,1,7,4,6,8),
  num.mark = c(4,2,4,1,3,2,2,2,1,3,1,1,4),
  seg.mean = c(5.3,2.3,1.2,1.1,1.4,2.2,3.3,0.5,3.3,4.3,4.3,6.3,7.3),
  Sample=c(rep("K",3),rep("L",5),rep("M",5)),
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
  checkEquals( segs2Rle( basic.segs[[1]], locData.rd ), basic.rle.df[[1]], checkNames=FALSE )
  checkEquals( segs2Rle( basic.segs[[2]], locData.rd ), basic.rle.df[[2]], checkNames=FALSE )
  checkEquals( segs2Rle( basic.segs[[3]], locData.rd ), basic.rle.df[[3]], checkNames=FALSE )
}

test_segs2RleDataFrame <- function() {
  checkEquals( segs2RleDataFrame( basic.segs, locData.rd ), basic.rle.df, checkNames=FALSE )
}

test_segTable <- function() {
  chr.ind = chrIndices(locData.rd)
  start = start(locData.rd)
  end = end(locData.rd)
  
  checkEquals( segTable( basic.rle.df[["K"]], locData.rd), basic.segs.after[["K"]], checkNames=FALSE )
  checkEquals( segTable( basic.rle.df[["L"]], locData.rd), basic.segs.after[["L"]], checkNames=FALSE )
  checkEquals( segTable( basic.rle.df[["M"]], locData.rd), basic.segs.after[["M"]], checkNames=FALSE )
  checkEquals( segTable( basic.rle.df[["M"]], chr.ind=chr.ind, start=start, end=end), basic.segs.after[["M"]], checkNames=FALSE, "segTable on Rle providing chr.ind, start, end" )
  checkEquals( segTable( basic.rle.df, locData.rd ), basic.segs.after, checkNames=FALSE )
  checkEquals( segTable( basic.rle.df, locData.rd, stack=TRUE ), stacked.basic.segs.after, checkNames=FALSE )
}


test_runCBS <- function() {
  sample.names = paste("a",1:2,sep="")
  probe.names =  paste("p",1:30,sep="")
  ds = matrix(c(c(rep(5,20),rep(3,10)),c(rep(2,10),rep(7,10),rep(9,10))),ncol=2,dimnames=list(probe.names,sample.names))
  ds.with.na = matrix(c(c(rep(5,9),NA,rep(5,10),rep(3,10)),c(rep(2,10),rep(7,10),rep(9,10))),ncol=2,dimnames=list(probe.names,sample.names))
  locs = RangedData(ranges=IRanges(start=c(1:20,1:10),width=1,names=probe.names),space=paste("chr",c(rep(1,20),rep(2,10)),sep=""))

  seg.rle.result = DataFrame( a1 = Rle(c(rep(5,20),rep(3,10))), a2 = Rle(c(rep(2,10),rep(7,10),rep(9,10))), row.names=probe.names )
  seg.list.result = list(
    a1 = data.frame( ID=rep("a1",2), chrom=factor(c("chr1","chr2")), loc.start=c(1,1), loc.end=c(20,10), num.mark=c(20,10), seg.mean=c(5,3), stringsAsFactors=FALSE),
    a2 = data.frame( ID=rep("a2",3), chrom=factor(c("chr1","chr1","chr2")), loc.start=c(1,11,1), loc.end=c(10,20,10), num.mark=c(10,10,10), seg.mean=c(2,7,9), stringsAsFactors=FALSE)
    )
  

  checkEquals( runCBS(ds,locs, n.cores=8), seg.rle.result, "Return DF of Rle")
  checkEquals( runCBS(ds.with.na,locs, n.cores=8), seg.rle.result, "Return DF of Rle with some NA in starting data")
  checkEquals( runCBS(ds,locs, n.cores=8, return.segs=TRUE), seg.list.result, "Return seg dfs")
  checkEquals( runCBS(seg.rle.result,locs, n.cores=8), seg.rle.result, "Return seg dfs starting from DF of Rle (like mbaf)")
  checkEquals( runCBS(ds,locs, n.cores=8,alpha=0.01), seg.rle.result, "Runs OK with alpha at 0.01 (requires loading of data from DNAcopy)")
}
