library(RUnit)
library(genoset)

### TO DO
# tidy up genome order tests
# GRanges tests for subsetting

test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]

test_creation <- function() {

  colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
  cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  locs=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  bad.locs=GRanges(ranges=IRanges(start=c(5,6,10:7,1:4),width=1,names=probe.names[c(5,6,10:7,1:4)]),seqnames=c(rep("chr3",2),rep("chrX",4),rep("chr1",4)))  

  tom = GenoSet( rowRanges=locs, cn=cn, colData=colData )

  gs.from.ad = GenoSet( rowRanges=locs, assayData=assayDataNew(storage.mode="environment",cn=cn), colData=colData )

  rle.genoset = GenoSet(
    rowRanges=locs,
    lrr=DataFrame(K=Rle(1:10),L=Rle(11:20),M=Rle(21:30),row.names=probe.names),
    baf=DataFrame(K=Rle(31:40),L=Rle(41:50),M=Rle(51:60),row.names=probe.names),
    colData=colData)
  
  misordered.genoset = GenoSet( rowRanges=locs, cn=cn[ rev(probe.names), ], foo=cn[ rev(probe.names),], colData=colData[rev(test.sample.names),] )

  bad.rowRanges.genoset = GenoSet( rowRanges=bad.locs, cn=cn, foo=cn, colData=colData )

  checkTrue(validObject(tom),"Regular GenoSet")
  checkTrue(validObject(gs.from.ad),"GenoSet with provided assayData")
  checkTrue(validObject(rle.genoset),"GenoSet with Rle data")
  checkTrue(validObject(misordered.genoset),"Starting with some sample name and feature name misordering")
  checkTrue( identical(misordered.genoset[,,"cn"],cn) && identical(misordered.genoset[,,"foo"],cn))
  checkIdentical( colData, colData(misordered.genoset), "Misordered colData gets fixed" )
  checkTrue(validObject(bad.rowRanges.genoset), "Can fix rowRanges not in strict genome order")
  checkIdentical( toGenomeOrder(rowRanges(bad.rowRanges.genoset),strict=TRUE), rowRanges(bad.rowRanges.genoset), "badly ordered rowRanges gets fixed" )
}

test_featureNames <- function() {
  ds = GenoSet(
    rowRanges=GRanges(IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))), 
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  gr = rowRanges(ds)
  checkEquals( featureNames(ds), probe.names, "Get featureNames from GenoSet" )
  checkEquals( featureNames(ds), featureNames(featureData(ds)), "GenoSet FeatureNames match featureData featureNames" )
  checkEquals( featureNames(ds), rownames(fData(ds)), "GenoSet FeatureNames match fData rownames" )
  checkEquals( featureNames(ds), featureNames(gr), "Get featureNames from GRanges" )
  checkEquals( featureNames(ds), rownames(ds), "featureNames and rownames are the same thing for a GenoSet.")
  checkEquals( featureNames(gr), rownames(gr), "featureNames and rownames are the same thing for a GRanges.")
  checkEquals( names(gr), rownames(gr), "names and rownames are the same thing for a GRanges.")
  new.fnames = paste("f",featureNames(ds),sep="")
  featureNames(ds) = new.fnames
  featureNames(gr) = new.fnames
  checkEquals( featureNames(ds), new.fnames, "Set featureNames in GenoSet")
  checkEquals( featureNames(gr), new.fnames, "Set featureNames in GRanges")
  new.fnames = paste("g",featureNames(ds),sep="")
  rownames(ds) = new.fnames
  rownames(gr) = new.fnames
  checkEquals( rownames(ds), new.fnames, "Set rownames in GenoSet")
  checkEquals( rownames(gr), new.fnames, "Set rownames in GRanges")
}

test_sampleNames <- function() {
  ds = GenoSet(
    rowRanges=GRanges(IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  checkIdentical( suppressWarnings(sampleNames(ds)), test.sample.names )
  checkIdentical( colnames(ds), test.sample.names )
  colnames(ds) = LETTERS[1:3]
  checkIdentical( colnames(ds), LETTERS[1:3] )
}

test_rowRanges <- function() {
  # With GRanges
  locs.gr = GRanges(ranges=IRanges(start=1:10, width=1, names=probe.names), seqnames=factor(c(rep("chr1", 4), rep("chr3", 2), rep("chrX", 4)), levels=c("chr1", "chr3", "chrX")))
  locs.gr.new = GRanges(ranges=IRanges(start=1:10, width=1, names=probe.names), seqnames=factor(c(rep("chr1", 4), rep("chr3", 2), rep("chrX", 4)), levels=c("chr1", "chr3", "chrX")))
  locs.gr.bad = locs.gr.new
  names(locs.gr.bad)[1] = "FOO"
  ds = GenoSet(
    rowRanges=locs.gr,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  ds.new = GenoSet(
    rowRanges=locs.gr.new,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  checkEquals(locs.gr,rowRanges(ds))
  checkEquals(ds,ds.new,check.attributes=FALSE)
  checkException( {rowRanges(ds) = locs.gr.bad},silent=TRUE )
  checkTrue({rowRanges(ds) = locs.gr; is.null(names(ds@rowRanges))}, "Setting rowRanges makes rowRanges GRanges names null")
  checkEquals(rownames(rowRanges(ds)), names(locs.gr), "However, getting rowRanges back out resets the GRanges names")
}

test_getters.and.setters <- function() {
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  colData = data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
  bob = GenoSet(
    rowRanges=GRanges(IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))), 
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    colData=colData
    )
  checkIdentical(c("baf", "lrr"), names(bob))
}

test_rd.gs.shared.api.and.getting.genome.info <- function() {
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  point.rowRanges = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  point.rowRanges.gr = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  point.bad.chr.order.rowRanges = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr5",4),rep("chrX",2),rep("chr3",4)))
  wide.rowRanges =  GRanges(ranges=IRanges(start=seq(1,30,by=3),width=3,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  gs = GenoSet(
    rowRanges=point.rowRanges,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  gr = as(point.rowRanges,"GRanges")

  checkEquals( start( point.rowRanges ), start( gs ) )
  checkEquals( width( point.rowRanges ), width( gs ) )
  checkEquals( end( point.rowRanges ), end( gs ) )
  checkEquals( chr(point.rowRanges), c(rep("chr1",4),rep("chr3",2),rep("chrX",4)) )
  checkEquals( chr(gr), c(rep("chr1",4),rep("chr3",2),rep("chrX",4)) )
  checkEquals( chr( point.rowRanges ), chr( gs ) )
  checkEquals( chr( point.rowRanges ), chr( gr ) )
  checkEquals( pos(point.rowRanges), 1L:10L )
  checkEquals( pos(wide.rowRanges), seq(from=2L, length=10, by=3L ) )
  checkEquals( pos( point.rowRanges ), pos( gs ) )
  checkEquals( pos( point.rowRanges ), pos( gr ) )
  checkEquals( chrNames( point.rowRanges ), c("chr1","chr3","chrX") )
  checkEquals( chrNames( point.rowRanges ), chrNames( gs ) )
  checkEquals( chrNames( gr[1:3,] ), c("chr1"), "chrNames on GRanges with empty levels should give just unique values" )
  point.rowRanges2 = point.rowRanges
  chrNames(point.rowRanges2) = sub("chr","",chrNames(point.rowRanges2))
  checkEquals( chrNames( point.rowRanges2 ), c("1","3","X") )
  gs2 = gs
  rowRanges(gs2) = point.rowRanges2
  checkEquals( chrNames( gs2 ), c("1","3","X") )
  checkEquals( chrNames( point.rowRanges ), c("chr1","chr3","chrX") )
  checkEquals( elementLengths( point.rowRanges ), elementLengths( gs ) )
  checkEquals( elementLengths( point.rowRanges ), elementLengths( point.rowRanges.gr ) )
  checkEquals( chrInfo( point.rowRanges ), chrInfo( gs ) )
  checkEquals( chrInfo( point.rowRanges ), chrInfo( gr ) )
  checkEquals( chrInfo( point.rowRanges ), matrix(c(1,5,11,4,10,20,0,4,10),ncol=3,dimnames=list(c("chr1","chr3","chrX"),c("start","stop","offset") ) ))
  checkEquals( chrIndices( point.rowRanges, "chr3"), c(5,6) )
  checkException( chrIndices( point.rowRanges, "chrFOO"), silent=TRUE )
  checkEquals( chrIndices( point.rowRanges ), chrIndices( gs ) )
  checkEquals( chrIndices( point.rowRanges ), matrix(c(1,5,7,4,6,10,0,4,6),ncol=3,dimnames=list(c("chr1","chr3","chrX"),c("first","last","offset") ) ))
  checkEquals( chrIndices( point.rowRanges ), chrIndices(point.rowRanges.gr) )
  checkEquals( chrIndices( point.rowRanges[1:6,] ), chrIndices(point.rowRanges.gr)[1:2,], "Empty levels ignored" )
  checkEquals( genoPos( point.rowRanges ), genoPos( gs ) )
  checkEquals( genoPos( point.rowRanges ), genoPos( gs ) )

  # Universe
  gr.uni = GRanges(IRanges(start=1:4,width=1),seqnames=c("chr1","chr2","chr3","chr4"))
  genome(gr.uni) = c("hg18","hg19","hg19","hg19")
  genome(gr.uni) = c("hg19")
  geno = rep("hg19", 3)
  rowRanges(gs) = gr
  genome(gs) = geno
  checkEquals(geno, genome(gs), "Get and set genome of GenoSet", checkNames=FALSE)
}

test_subset <- function() {
  test.gr = GRanges(ranges=IRanges(start=8:14,width=1,names=letters[8:14]),seqnames=rep("chrX",7))
  test.pdata = data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])),stringsAsFactors=FALSE)
  test.phenodata = phenoData=new("AnnotatedDataFrame",test.pdata)

  test.ds = new("GenoSet",
    rowRanges=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    phenoData=test.phenodata
    )
  
  expected.ds = new("GenoSet",
    rowRanges=GRanges(ranges=IRanges(start=8:10,width=1,names=probe.names[8:10]),seqnames="chrX"),
    lrr=matrix(c(8:10,18:20,28:30),nrow=3,ncol=3,dimnames=list(probe.names[8:10],test.sample.names)),
    baf=matrix(c(38:40,48:50,58:60),nrow=3,ncol=3,dimnames=list(probe.names[8:10],test.sample.names)),
    phenoData=test.phenodata
    )
  
  chr3.ds = new("GenoSet",
    rowRanges=GRanges(ranges=IRanges(start=5:6,width=1,names=probe.names[5:6]),seqnames="chr3"),
    lrr=matrix(c(5:6,15:16,25:26),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
    baf=matrix(c(35:36,45:46,55:56),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
    phenoData=test.phenodata
    )

  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  
  ds = GenoSet(
    rowRanges=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    colData=test.pdata,
    annotation="SNP6"
    )

  subset.rows.ds = GenoSet(
    rowRanges=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))[2:3,,drop=TRUE],
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))[2:3,],
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))[2:3,],
    colData=test.pdata,
    annotation="SNP6"
    )
  
  subset.cols.ds = GenoSet(
    rowRanges=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    lrr=matrix(11:30,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
    baf=matrix(41:60,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
    colData=test.pdata[2:3,],
    annotation="SNP6"
    )

  gene.gr = GRanges(ranges=IRanges(start=2:3,width=1),seqnames=c("chr1","chr1"))
  
  # Subsetting whole object
  checkEquals( ds[ ,2:3], subset.cols.ds, check.attributes=TRUE)
  checkEquals( ds[ 2:3, ], subset.rows.ds, check.attributes=FALSE)
  checkEquals( ds[ 2:3, , "baf" ], subset.rows.ds[, , "baf"], check.attributes=TRUE)
  ds.two.rows = ds[ 2:3, ]
  checkEquals( featureNames(assayData(ds.two.rows)), featureNames(rowRanges(ds.two.rows)), "featureNames from rowRanges and assayData should be the same when rows subset.")
  checkEquals( ds[ gene.gr, ], subset.rows.ds, check.attributes=FALSE )
  
  # Subsetting assayData / extracting
  checkEquals( ds[ 5, 3, "baf"], assayDataElement(ds,"baf")[5,3])
  checkEquals( ds[ 5, 3, 1], assayDataElement(ds,"baf")[5,3])
  checkEquals( ds[ , , "lrr"], assayDataElement(ds,"lrr"), "Extract whole matrix" )
  checkException( ds[ , , "foo"], "Fail to extract assayDataElement with bad character k", silent=TRUE)
  checkException( ds[ , , 8], "Fail to extract assayDataElement with bad integer k", silent=TRUE)
  
  # Test subsetting by location
  checkEquals( test.ds[test.gr,], expected.ds, check.attributes=FALSE )
  checkEquals( test.ds[as(test.gr,"GRanges"),], expected.ds, check.attributes=FALSE )
  checkEquals( test.ds[8:10,], expected.ds, check.attributes=FALSE)
  checkEquals( test.ds[ chrIndices(test.ds,"chr3"), ], chr3.ds , check.attributes=FALSE )

  # Replace
  ds = GenoSet(
    rowRanges=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )

  ds[,,"baf"] = ds[,,"lrr"]
  checkEquals(ds[,,"baf"],ds[,,"lrr"],"Replace whole element")
  bad.names.lrr = ds[,,"lrr"]
  rownames(bad.names.lrr)[1] = "FOO"
  colnames(bad.names.lrr)[1] = "FOO"
  checkException({ds[,,"baf"] = bad.names.lrr}, "Incoming ad element must have dimnames that matches genoset.", silent=TRUE)
  lrr.mat = ds[,,"lrr"]
  lrr.mat[1:2,1:2] = 5
  ds[1:2,1:2,"lrr"] = 5
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace partial matrix with integer indices")
  lrr.mat[6:8,2] = 3
  ds[rowRanges(ds)[6:8,],2,"lrr"] = 3
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace partial matrix with GRanges subsetting of rows")
  ds[,3,"lrr"] = 3
  lrr.mat[,3] = 3
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace column")
  ds[, , "FOO"] = ds[, , "lrr"]
  checkEquals(ds[, , "FOO"], ds[, , "lrr"], "Adding a whole new matrix is OK.")
  checkException({ds[1, 1, "foo"] = 5}, "Fail to replace with bad character assayDataElement index k", silent=TRUE)
  checkException({ds[1, 1, 8] = 5}, "Fail to replace with bad integer assayDataElement index k", silent=TRUE)
}

test_subset_w_granges <- function() {
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  test.gr = GRanges(ranges=IRanges(start=8:14,width=1),names=letters[8:14],seqnames=rep("chrX",7))
  lrr = matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  baf = matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  colData = data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
  locs = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  test.ds = GenoSet(
    rowRanges=locs,
    lrr=lrr,
    baf=baf,
    colData=colData
    )
  
  expected.ds = GenoSet(
    rowRanges=locs[8:10,],
    lrr=lrr[8:10,],
    baf=baf[8:10,],
    colData=colData
    )
  
  chr3.ds = GenoSet(
    rowRanges=locs[5:6,],
    lrr=matrix(c(5:6,15:16,25:26),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
    baf=matrix(c(35:36,45:46,55:56),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
    colData=colData
    )
  
  ds = GenoSet(
    rowRanges=locs,
    lrr=lrr,
    baf=baf,
    colData=colData
    )

  subset.rows.ds = GenoSet(
    rowRanges=locs[2:3,],
    lrr=lrr[2:3,],
    baf=baf[2:3,],
    colData=colData
    )
  
  subset.cols.ds = GenoSet(
    rowRanges=locs,
    lrr=matrix(11:30,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
    baf=matrix(41:60,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
    colData=colData[2:3,]
    )

  gene.gr = GRanges(ranges=IRanges(start=2:3,width=1),seqnames=c("chr1","chr1"))
  
  # Subsetting whole object
  checkEquals( ds[ ,2:3], subset.cols.ds, check.attributes=FALSE)
  checkEquals( ds[ 2:3, ], subset.rows.ds, check.attributes=FALSE)
  checkEquals( ds[ gene.gr, ], subset.rows.ds, check.attributes=FALSE)
  
  # Subsetting assayData / extracting
  checkEquals( ds[ 5, 3, "baf"], assay(ds,"baf")[5,3])
  checkEquals( ds[ 5, 3, 1], assay(ds,"baf")[5,3])
  checkEquals( ds[ , , "lrr"], assay(ds,"lrr"), "Extract whole matrix" )
  
  # Test subsetting by location
  checkEquals( test.ds[test.gr,], expected.ds, check.attributes=FALSE)
  checkEquals( test.ds[8:10,], expected.ds, check.attributes=FALSE)
  checkEquals( test.ds[ chrIndices(test.ds,"chr3"), ], chr3.ds, check.attributes=FALSE)

  # Replace
  ds = GenoSet(
    rowRanges=locs,
    lrr=lrr,
    baf=baf,
    colData=colData
    )

  ds[,,"baf"] = ds[,,"lrr"]
  checkEquals(ds[,,"baf"],ds[,,"lrr"],"Replace whole element")
  bad.names.lrr = ds[,,"lrr"]
  rownames(bad.names.lrr)[1] = "FOO"
  colnames(bad.names.lrr)[1] = "FOO"
  checkException({ds[,,"baf"] = bad.names.lrr}, "Incoming ad element must have dimnames that matches genoset.",silent=TRUE)
  lrr.mat = ds[,,"lrr"]
  lrr.mat[1:2,1:2] = 5
  ds[1:2,1:2,"lrr"] = 5
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace partial matrix with integer indices")
  lrr.mat[6:8,2] = 3
  ds[rowRanges(ds)[6:8,],2,"lrr"] = 3
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace partial matrix with RangedData subsetting of rows")
  ds[,3,"lrr"] = 3
  lrr.mat[,3] = 3
  checkEquals(lrr.mat,ds[,,"lrr"],"Replace column")
}

test_genomeOrder <- function() {
  chr.names = c(rep("chr1",3),rep("chr2",3),rep("chr10",4))

  ok.locs = GRanges( ranges = IRanges(start=1:10,width=1,names=paste("p",1:10,sep="")), seqnames=factor(chr.names,levels=c("chr1","chr2","chr10")))
  checkTrue( isGenomeOrder(ok.locs), "Good locs" )

  bad.locs = GRanges( ranges = IRanges(start=c(2,3,1,4,6,5,10:7),width=1,names=paste("p",c(2,3,1,4,6,5,10:7),sep="")), seqnames=factor(chr.names,levels=c("chr1","chr2","chr10")))
  bad.locs.bad.chr = GRanges( ranges = IRanges(start=c(2,3,1,4,6,5,10:7),width=1,names=paste("p",c(2,3,1,4,6,5,10:7),sep="")), seqnames=factor(chr.names,levels=c("chr2","chr1","chr10")))
  checkTrue( ! isGenomeOrder(bad.locs, strict=TRUE), "Bad within chr, OK chr levels, fail")

  good.ds = GenoSet(
    rowRanges=ok.locs,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(rownames(ok.locs),test.sample.names)),
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  bad.ds = good.ds[ c(6,5,4,3,2,1,10:7),]
  bad.ds.bad.chrs = GenoSet(
    rowRanges=bad.locs.bad.chr,
    cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(rownames(ok.locs),test.sample.names))[c(4,6,5,2,3,1,10:7),],
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  checkTrue(isGenomeOrder(good.ds))
  checkTrue(!isGenomeOrder(bad.ds))
  checkEquals( good.ds, toGenomeOrder(bad.ds,strict=TRUE), check.attributes=FALSE, "GenoSet disordered within chrs" )
  checkEquals( good.ds, toGenomeOrder(bad.ds.bad.chrs,strict=TRUE), check.attributes=FALSE, "GenoSet disordered within chrs, disordered chrs" )

  gr1 = GRanges(ranges=IRanges(start=c(9,1,5,4,6,2),width=1,names=LETTERS[c(9,1,5,4,6,2)]),seqnames=Rle(factor(c("A","B","C","C","B","A"),levels=c("A","C","B"))))
  gr2 = GRanges(ranges=IRanges(start=c(2,9,4,5,1,6),width=1,names=LETTERS[c(2,9,4,5,1,6)]),seqnames=Rle(factor(c("A","A","C","C","B","B"),levels=c("A","C","B"))))
  gr3 = GRanges(ranges=IRanges(start=c(2,9,1,6,4,5),width=1,names=LETTERS[c(2,9,1,6,4,5)]),seqnames=Rle(factor(c("A","A","B","B","C","C"),levels=c("A","B","C"))))
  checkIdentical(toGenomeOrder(gr1,strict=FALSE),gr2,"GRanges with mis-ordered chromosomes, without strict")
  checkIdentical(toGenomeOrder(gr1,strict=TRUE),gr3,"GRanges with mis-ordered chromosomes, with strict")
  checkTrue(isGenomeOrder(gr2,strict=FALSE))
  checkTrue(isGenomeOrder(gr3,strict=TRUE))
  checkTrue(!isGenomeOrder(gr2,strict=TRUE))
  checkTrue(!isGenomeOrder(gr1,strict=TRUE), "Not in blocks by chromsome, strict")
  checkTrue(!isGenomeOrder(gr1,strict=FALSE), "Not in blocks by chromsome, strict")
}
