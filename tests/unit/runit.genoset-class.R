test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]

test.creation <- function() {

  bob = BAFSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  locData.rd = RangedData(ranges=IRanges(start=c(1,4,3,2,5:10),width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18")
  
  ted = BAFSet(
    locData=locData.rd,
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  joe = CNSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  tom = GenoSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  # Test making a BAFSet with Rle-based data
  rle.bafset = BAFSet(
    locData=locData.rd,
    lrr=DataFrame(K=Rle(1:10),L=Rle(11:20),M=Rle(21:30),row.names=probe.names),
    baf=DataFrame(K=Rle(31:40),L=Rle(41:50),M=Rle(51:60),row.names=probe.names),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  # Test making a CNSet with Rle-based data
  rle.cnset = CNSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    cn=DataFrame(K=Rle(1:10),L=Rle(11:20),M=Rle(21:30),row.names=probe.names),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  
  checkTrue(validObject(bob),"Regular BAFSet")
  checkTrue(validObject(ted),"BAFSet with out of genome order locData")
  checkTrue(validObject(joe),"CNSet with out of genome order locData")
  checkTrue(validObject(tom),"GenoSet with out of genome order locData")
  checkTrue(validObject(rle.bafset),"BAFSet with Rle data")
  checkTrue(validObject(rle.cnset),"CNSet with Rle-based data")
}

test.sampleNames <- function() {
  ds = CNSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  bad.sampleNames = c("K-fed","&FOO","&FOO")
  checkEquals( sampleNames(ds), test.sample.names )
  sampleNames(ds) = bad.sampleNames
  checkEquals( sampleNames(ds), c("K.fed","X.FOO","X.FOO.1") )
}

test.featureNames <- function() {
  ds = CNSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  bad.featureNames = c("a.","b,,","-foo",letters[4:9],"b,,")
  checkEquals( featureNames(ds), probe.names )
  featureNames(ds) = bad.featureNames
  checkEquals( featureNames(ds), c("a.","b..","X.foo",letters[4:9],"b...1"))
}

test.locData <- function() {
  ld = RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=factor(c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),levels=c("chr1","chr3","chrX")),universe="hg18")
  ds = CNSet(
    locData=ld,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  checkEquals(ld,locData(ds))
  ld.new = RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=factor(c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),levels=c("chr1","chr3","chrX")),universe="hg18")
  locData(ds) = ld.new
  ds.new = CNSet(
    locData=ld.new,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  checkEquals(ds,ds.new,check.attributes=FALSE)
  ld.bad = ld.new
  rownames(ld.bad)[1] = "FOO"
  checkException( eval(parse(text="locData(ds) = ld.bad")) )
}

test.rd.gs.shared.api.and.getting.genome.info <- function() {
  point.locData = RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg19")
  point.bad.chr.order.locData = RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr5",4),rep("chrX",2),rep("chr3",4)),universe="hg19")
  wide.locData =  RangedData(ranges=IRanges(start=seq(1,30,by=3),width=3,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg19")
  gs = GenoSet(
    locData=point.locData,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  checkEquals( chr(point.locData), c(rep("chr1",4),rep("chr3",2),rep("chrX",4)) )
  checkEquals( chr( point.locData ), chr( gs ) )
  checkEquals( pos(point.locData), 1L:10L )
  checkEquals( pos(wide.locData), seq(from=2L, length=10, by=3L ) )
  checkEquals( pos( point.locData ), pos( gs ) )
  checkEquals( uniqueChrs( point.locData ), c("chr1","chr3","chrX") )
  checkEquals( uniqueChrs( point.locData ), uniqueChrs( gs ) )
  checkEquals( names( point.locData ), c("chr1","chr3","chrX") )
  checkEquals( names( point.locData ), names( gs ) )
  checkEquals( ranges( point.locData ), ranges( gs ) )
  checkEquals( elementLengths( point.locData ), elementLengths( gs ) )
  checkEquals( orderedChrs( point.bad.chr.order.locData ), c("chr3","chr5","chrX") )
  checkEquals( orderedChrs( point.locData ), orderedChrs( gs ) )
  checkEquals( chrInfo( point.locData ), chrInfo( gs ) )
  checkEquals( chrInfo( point.locData ), matrix(c(1,5,11,4,10,20,0,4,10),ncol=3,dimnames=list(c("chr1","chr3","chrX"),c("start","stop","offset") ) ))
  checkEquals( chrIndices( point.locData, "chr3"), c(5,6) )
  checkException( chrIndices( point.locData, "chrFOO") )
  checkEquals( chrIndices( point.locData ), chrIndices( gs ) )
  checkEquals( chrIndices( point.locData ), matrix(c(1,5,7,4,6,10,0,4,6),ncol=3,dimnames=list(c("chr1","chr3","chrX"),c("first","last","offset") ) ))
  checkEquals( genoPos( point.locData ), genoPos( gs ) )
  checkEquals( genoPos( point.locData ), genoPos( gs ) )
}

test.subset <- function() {
  
  test.rd = RangedData(ranges=IRanges(start=8:14,width=1),names=letters[8:14],space=rep("chrX",7))
    
  test.ds = new("BAFSet",
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    phenoData=new("AnnotatedDataFrame",data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))))
    )
  
  expected.ds = new("BAFSet",
    locData=RangedData(ranges=IRanges(start=8:10,width=1,names=probe.names[8:10]),space="chrX",universe="hg18"),
    lrr=matrix(c(8:10,18:20,28:30),nrow=3,ncol=3,dimnames=list(probe.names[8:10],test.sample.names)),
    baf=matrix(c(38:40,48:50,58:60),nrow=3,ncol=3,dimnames=list(probe.names[8:10],test.sample.names)),
    phenoData=new("AnnotatedDataFrame",data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))))
    )
  
  chr3.ds = new("BAFSet",
    locData=RangedData(ranges=IRanges(start=5:6,width=1,names=probe.names[5:6]),space="chr3",universe="hg18"),
    lrr=matrix(c(5:6,15:16,25:26),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
    baf=matrix(c(35:36,45:46,55:56),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
    phenoData=new("AnnotatedDataFrame",data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))))
    )

  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  
  ds = BAFSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  subset.rows.ds = BAFSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18")[2:3,,drop=TRUE],
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))[2:3,],
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))[2:3,],
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  
  subset.cols.ds = BAFSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    lrr=matrix(11:30,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
    baf=matrix(41:60,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))[2:3,]),
    annotation="SNP6"
    )

  gene.rd = RangedData(ranges=IRanges(start=2:3,width=1),space=c("chr1","chr1"),universe="hg18")

  bm.ds = convertToBigMatrix(ds,path=tempdir())
  
  # Subsetting whole object
  checkEquals( ds[ ,2:3], subset.cols.ds, check.attributes=FALSE)
  checkEquals( ds[ 2:3, ], subset.rows.ds, check.attributes=FALSE)
  checkEquals( ds[ gene.rd, ], subset.rows.ds, check.attributes=FALSE)
  
  # Subsetting assayData / extracting
  checkEquals( ds[ 4:6, 1:2, "lrr"], lrr(ds)[4:6,1:2])
  checkEquals( ds[ 5, 3, "baf"], baf(ds)[5,3])
  checkEquals( ds[ 5, 3, "baf"], assayDataElement(ds,"baf")[5,3])
  checkEquals( ds[ 5, 3, 1], assayDataElement(ds,"baf")[5,3])
  checkEquals( ds[ gene.rd, 1:2,"lrr" ], lrr(ds)[2:3,1:2] )
  checkEquals( ds[ , , "lrr"], assayDataElement(ds,"lrr"), "Extract whole matrix" )
  checkEquals( bm.ds[ 1:3, 1:3, "lrr"], assayDataElement(bm.ds,"lrr")[1:3,1:3], "Extract part of big.matrix" )
  checkIdentical( bm.ds[ , , "lrr"], assayDataElement(bm.ds,"lrr"), "Extract whole big.matrix" )
  
  # Test subsetting by location
  checkEquals( test.ds[test.rd,], expected.ds, checkNames=FALSE )
  checkEquals( test.ds[ranges(test.rd),], expected.ds, checkNames=FALSE )
  checkEquals( test.ds[8:10,], expected.ds, checkNames=FALSE )
  checkEquals( test.ds[ chrIndices(test.ds,"chr3"), ], chr3.ds , checkNames=FALSE)
}

test.gcCorrect <- function() {

  input.vector = c(rep(0.05,50),rep(0.08,50))
  gc = input.vector
  output.vector = rep(0,100)
  checkEquals( gcCorrect(input.vector, gc ), output.vector )

  input.matrix = matrix(c(input.vector,input.vector),ncol=2)
  output.matrix = matrix(c(output.vector,output.vector),ncol=2)
  checkEquals( gcCorrect(input.matrix, gc ), output.matrix )
  
  gc.w.na = gc
  is.na(gc.w.na) = c(25,75)
  output.matrix.w.na = output.matrix
  output.matrix.w.na[ c(25,75), ] = NA
  checkEquals( gcCorrect(input.matrix, gc.w.na ), output.matrix.w.na )

}

test.genomeOrder <- function() {
  chr.names = c(rep("chr1",3),rep("chr2",3),rep("chr10",4))

  ok.locs = RangedData( ranges = IRanges(start=1:10,width=1,names=paste("p",1:10,sep="")), space=factor(chr.names,levels=c("chr1","chr2","chr10")))
  checkTrue( isGenomeOrder(ok.locs), "Good locs" )

  ok.locs.weak = RangedData( ranges = IRanges(start=1:10,width=1,names=paste("p",1:10,sep="")), space=factor(chr.names,levels=c("chr1","chr10","chr2")))
  checkTrue( isGenomeOrder(ok.locs.weak, strict=FALSE), "Good locs with disordered chrs OK" )
  checkTrue( ! isGenomeOrder(ok.locs.weak, strict=TRUE), "Good locs with disordered chrs not strict" )

  bad.locs = RangedData( ranges = IRanges(start=c(2,3,1,4,6,5,10:7),width=1,names=paste("p",c(2,3,1,4,6,5,10:7),sep="")), space=factor(chr.names,levels=c("chr1","chr2","chr10")))
  bad.locs.bad.chr = RangedData( ranges = IRanges(start=c(2,3,1,4,6,5,10:7),width=1,names=paste("p",c(2,3,1,4,6,5,10:7),sep="")), space=factor(chr.names,levels=c("chr2","chr1","chr10")))
  checkTrue( ! isGenomeOrder(bad.locs, strict=TRUE), "Bad within chr, OK chr levels, fail")

  checkEquals( ok.locs, toGenomeOrder(ok.locs,strict=TRUE), "Perfect locs pass")
  checkEquals( ok.locs.weak, toGenomeOrder(ok.locs.weak,strict=FALSE), "locs with disordered chr block pass with strict as FALSE")
  checkEquals( ok.locs, toGenomeOrder(ok.locs.weak,strict=TRUE), "locs disordered chrs, but ok within chrs passes with strict as TRUE")
  checkEquals( ok.locs, toGenomeOrder(bad.locs,strict=TRUE), "locs ok chrs, but disordered within")
  checkEquals( ok.locs, toGenomeOrder(bad.locs.bad.chr,strict=TRUE), "locs with disordered chrs and within chrs")

  good.ds = CNSet(
    locData=ok.locs,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(rownames(ok.locs),test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6", universe="hg19"
    )
  bad.ds = good.ds[ c(6,5,4,3,2,1,10:7),]
  bad.ds.bad.chrs = CNSet(
    locData=bad.locs.bad.chr,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(rownames(ok.locs),test.sample.names))[c(4,6,5,2,3,1,10:7),],
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6", universe="hg19"
    )
  checkTrue(isGenomeOrder(good.ds))
  checkTrue(!isGenomeOrder(bad.ds))
  checkEquals( good.ds, toGenomeOrder(bad.ds,strict=TRUE), check.attributes=FALSE, "CNSet disordered within chrs" )
  checkEquals( good.ds, toGenomeOrder(bad.ds.bad.chrs,strict=TRUE), check.attributes=FALSE, "CNSet disordered within chrs, disordered chrs" )
}

test.loadGC <- function() {

}
