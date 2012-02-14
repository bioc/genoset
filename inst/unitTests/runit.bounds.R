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
}
