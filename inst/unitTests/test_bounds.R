# Tests for functions utilizing boundingIndices

test_boundingIndices <- function() {

  # Test with exact matches
  gene.starts = seq( 0, 42, 2)
  gene.stops = gene.starts + 2
  probes = 1:40
  bounds = matrix( c(c(seq(0,40,2),40),c(seq(2,40,2),41,41)), ncol=2)
  valid.bounds = matrix( c(c(1,seq(2,40,2),40),c(seq(2,40,2),40,40)), ncol=2)

  checkEquals( boundingIndices(gene.starts, gene.stops, probes, valid.indices=FALSE), bounds)
  checkEquals( boundingIndices(gene.starts, gene.stops, probes, valid.indices=TRUE), valid.bounds)
  checkEquals( boundingIndices2(gene.starts, gene.stops, probes), valid.bounds)

  # Test random order input with some exact matches
  gene.starts = c(6,2,9,1,14,7,50)
  gene.stops = gene.starts + 2
  probes = seq(3,39,3)
  bounds = matrix( c(c(2,0,3,0,4,2,13),c(3,2,4,1,6,3,14)), ncol=2)
  valid.bounds = matrix( c(c(2,1,3,1,4,2,13),c(3,2,4,1,6,3,13)), ncol=2)

  checkEquals( boundingIndices(gene.starts, gene.stops, probes, valid.indices=FALSE), bounds)
  checkEquals( boundingIndices(gene.starts, gene.stops, probes, valid.indices=TRUE), valid.bounds)
  checkEquals( boundingIndices2(gene.starts, gene.stops, probes), valid.bounds)

  # Test with some not matching exactly
  gene.starts = seq(1,16,4)
  gene.stops = seq(1,16,4) + 1
  probes = seq(2,14,3)
  valid.bounds = matrix(c(1,2,3,4,1,3,4,5),ncol=2)
  bounds = matrix(c(0,2,3,4,1,3,4,5),ncol=2)
  
  checkEquals( boundingIndices(gene.starts, gene.stops, probes, valid.indices=TRUE), valid.bounds)
  checkEquals( boundingIndices(gene.starts, gene.stops, probes, valid.indices=FALSE), bounds)
  checkEquals( boundingIndices(gene.starts, gene.stops, probes, valid.indices=FALSE, offset=2), bounds + 2)
  checkEquals( boundingIndices2(gene.starts, gene.stops, probes), valid.bounds)
  
}

test_rangeSampleMeans <- function() {
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
 
  subject = CNSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  query.rd = RangedData( ranges=IRanges(start=c(2,3,7,8),width=2,names=c("joe","bob","fred","tom")), space=factor(c("chr1","chr1","chrX","chrX"),levels=c("chr1","chrX")),universe="hg18")

  means = matrix(c(32,42,52,33,43,53,37,47,57,38,48,58)+0.5,ncol=nrow(query.rd),nrow=ncol(subject),dimnames=list(sampleNames(subject),rownames(query.rd)))
  means = t(means)
  checkEquals( rangeSampleMeans( query.rd, subject, "cn" ), means)

  bigmat.dir = file.path(tempdir(),"bigmat")
  bm.ds = convertToBigMatrix(subject,path=bigmat.dir)
  checkEquals( rangeSampleMeans( query.rd, bm.ds, "cn" ), means)
  rm.results = try(unlink(bigmat.dir,recursive=TRUE),silent=TRUE)
  checkTrue( !inherits(rm.results,"try-error") )
  
  rle.cnset = CNSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    cn=DataFrame(K=Rle(1:10),L=Rle(11:20),M=Rle(21:30),row.names=probe.names),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  rle.means = matrix(c(2.5,3.5,7.5,8.5,12.5,13.5,17.5,18.5,22.5,23.5,27.5,28.5), nrow=nrow(query.rd), ncol=ncol(rle.cnset), dimnames=list(rownames(query.rd),sampleNames(rle.cnset)))
  checkEquals( rangeSampleMeans( query.rd, rle.cnset, "cn" ), rle.means, "DataFrame of Rle")
}

test_rangeColMeans <- function() {
  bounds = matrix(c(2,3,3,5,7,8,9,10),ncol=2,byrow=TRUE)
  x = matrix(31:60,nrow=10,ncol=3)
  means = matrix(c(32.5,34,37.5,39.5,42.5,44,47.5,49.5,52.5,54,57.5,59.5),nrow=nrow(bounds),ncol=ncol(x))
  checkEquals( rangeColMeans( bounds, x), means, "Matrix without dimnames")
  checkEquals( rangeColMeans( bounds, x[,1]), means[,1], "Vector without dimnames")
  
  rownames(x) = letters[1:nrow(x)]
  colnames(x) = letters[1:ncol(x)]
  rownames(bounds) = LETTERS[1:nrow(bounds)]
  rownames(means) = rownames(bounds)
  colnames(means) = colnames(x)
  checkEquals( rangeColMeans( bounds, x), means, "Matrix with dimnames")
  checkEquals( rangeColMeans( bounds, x[,1]), means[,1], "Vector without dimnames")

  na.cells = matrix(c(3,1,4,2,8,1,8,3,2,3,3,3),ncol=2,byrow=TRUE)
  x.w.na = x
  x.w.na[ na.cells ] = NA
  x.w.na[ 8,3 ] = NaN
  means.w.na = matrix(c(32,34.5,37,39.5, 42.5,44,47.5,49.5, NA,54.5,57,59.5),nrow=nrow(bounds),ncol=ncol(x),dimnames=dimnames(means))
  checkEquals( rangeColMeans( bounds, x.w.na), means.w.na, "Matrix with dimnames and NAs")
}

test_boundingIndicesByChr <- function() {
  subject= RangedData(ranges=IRanges(start=c(seq(from=10,to=40,by=10),seq(from=110,to=140,by=10),seq(from=1110,to=1140,by=10)),width=2,names=as.character(1:12)),space=c(rep("1",4),rep("2",4),rep("3",4)))
  query = RangedData(ranges=IRanges(start=c(2,9,39,50,102,109,139,150,1102,1109,1139,1150),width=2,names=as.character(1:12)),space=c(rep("1",4),rep("2",4),rep("3",4)))
  res = matrix(as.integer(c(1,1, 1,1, 3,4, 4,4, 5,5, 5,5, 7,8, 8,8, 9,9, 9,9, 11,12, 12,12)),byrow=TRUE,ncol=2,dimnames=list(rownames(query),c("left","right")))
  checkIdentical(res, boundingIndicesByChr(query,subject))

  subject2= RangedData(ranges=IRanges(start=c(seq(from=10,to=40,by=10),seq(from=110,to=140,by=10),seq(from=1110,to=1140,by=10)),width=2,names=as.character(1:12)),space=c(rep("1",4),rep("2",4),rep("5",4)))
  query2 = RangedData(ranges=IRanges(start=c(2,9,39,50,102,109,139,150,1102,1109,1139,1150),width=2,names=as.character(1:12)),space=c(rep("1",4),rep("3",4),rep("5",4)))
  res2 = res[c(1:4,9:12),]
  checkIdentical(res2, boundingIndicesByChr(query2,subject2))

}
