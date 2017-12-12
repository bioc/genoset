library(testthat)
library(genoset)

### TO DO
# tidy up genome order tests
# GRanges tests for subsetting

test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]

test_that("We can make genosets", {

  colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
  cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  locs=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  bad.locs=GRanges(ranges=IRanges(start=c(5,6,10:7,1:4),width=1,names=probe.names[c(5,6,10:7,1:4)]),seqnames=c(rep("chr3",2),rep("chrX",4),rep("chr1",4)))

  tom = GenoSet( rowRanges=locs, assays=list(cn=cn), colData=colData )

  rle.genoset = GenoSet(
    rowRanges=locs,
    assays=list(
        lrr=RleDataFrame(K=Rle(1:10),L=Rle(11:20),M=Rle(21:30),row.names=probe.names),
        baf=RleDataFrame(K=Rle(31:40),L=Rle(41:50),M=Rle(51:60),row.names=probe.names)
        ),
    colData=colData)

    expect_error(
        GenoSet( rowRanges=locs, assays=list(cn=cn[ rev(probe.names), ], foo=cn[ rev(probe.names),]),
                colData=colData[rev(test.sample.names),]
                ),
        silent=TRUE
        )

  expect_true(validObject(tom),"Regular GenoSet")
  expect_true(validObject(rle.genoset),"GenoSet with Rle data")
})

test_that("Our shared GRanges genoset API works", {
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  point.rowRanges = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  point.rowRanges.gr = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  point.bad.chr.order.rowRanges = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr5",4),rep("chrX",2),rep("chr3",4)))
  wide.rowRanges =  GRanges(ranges=IRanges(start=seq(1,30,by=3),width=3,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  gs = GenoSet(
    rowRanges=point.rowRanges,
    assays=list(
        cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
        ),
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
    )
  gr = as(point.rowRanges,"GRanges")

  expect_equal( start( point.rowRanges ), start( gs ) )
  expect_equal( width( point.rowRanges ), width( gs ) )
  expect_equal( end( point.rowRanges ), end( gs ) )
  expect_equal( chr(point.rowRanges), c(rep("chr1",4),rep("chr3",2),rep("chrX",4)) )
  expect_equal( chr(gr), c(rep("chr1",4),rep("chr3",2),rep("chrX",4)) )
  expect_equal( chr( point.rowRanges ), chr( gs ) )
  expect_equal( chr( point.rowRanges ), chr( gr ) )
  expect_equal( pos(point.rowRanges), 1L:10L )
  expect_equal( pos(wide.rowRanges), seq(from=2L, length=10, by=3L ) )
  expect_equal( pos( point.rowRanges ), pos( gs ) )
  expect_equal( pos( point.rowRanges ), pos( gr ) )
  expect_equal( chrNames( point.rowRanges ), c("chr1","chr3","chrX") )
  expect_equal( chrNames( point.rowRanges ), chrNames( gs ) )
  expect_equal( chrNames( gr[1:3,] ), c("chr1"), label = "chrNames on GRanges with empty levels should give just unique values" )
  point.rowRanges2 = point.rowRanges
  chrNames(point.rowRanges2) = sub("chr","",chrNames(point.rowRanges2))
  expect_equal( chrNames( point.rowRanges2 ), c("1","3","X") )
  gs2 = gs
  rowRanges(gs2) = point.rowRanges2
  expect_equal( chrNames( point.rowRanges ), c("chr1","chr3","chrX") )
  expect_equal( lengths( point.rowRanges ), lengths( gs ) )
  expect_equal( lengths( point.rowRanges ), lengths( point.rowRanges.gr ) )
  expect_equal( chrInfo( point.rowRanges ), chrInfo( gs ) )
  expect_equal( chrInfo( point.rowRanges ), chrInfo( gr ) )
  expect_equal( chrInfo( point.rowRanges ), matrix(c(1,5,11,4,10,20,0,4,10),ncol=3,dimnames=list(c("chr1","chr3","chrX"),c("start","stop","offset") ) ))
  expect_equal( chrIndices( point.rowRanges, "chr3"), c(5,6) )
  expect_error( chrIndices( point.rowRanges, "chrFOO"), silent=TRUE )
  expect_equal( chrIndices( point.rowRanges ), chrIndices( gs ) )
  expect_equal( chrIndices( point.rowRanges ), matrix(c(1,5,7,4,6,10,0,4,6),ncol=3,dimnames=list(c("chr1","chr3","chrX"),c("first","last","offset") ) ))
  expect_equal( chrIndices( point.rowRanges ), chrIndices(point.rowRanges.gr) )
  expect_equal( chrIndices( point.rowRanges[1:6,] ), chrIndices(point.rowRanges.gr)[1:2,], label = "Empty levels ignored" )
  expect_equal( genoPos( point.rowRanges ), genoPos( gs ) )
  expect_equal( genoPos( point.rowRanges ), genoPos( gs ) )
})

test_that("We can subset genoset", {
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  test.gr = GRanges(ranges=IRanges(start=8:14,width=1),names=letters[8:14],seqnames=rep("chrX",7))
  lrr = matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  baf = matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names))
  colData = data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
  locs = GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
  test.ds = GenoSet(
    rowRanges=locs,
    assays=list(
        lrr=lrr,
        baf=baf
        ),
    colData=colData
    )

  expected.ds = GenoSet(
    rowRanges=locs[8:10,],
    assays=list(
        lrr=lrr[8:10,],
        baf=baf[8:10,]
        ),
    colData=colData
    )

  chr3.ds = GenoSet(
    rowRanges=locs[5:6,],
    assays=list(
        lrr=matrix(c(5:6,15:16,25:26),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names)),
        baf=matrix(c(35:36,45:46,55:56),nrow=2,ncol=3,dimnames=list(probe.names[5:6],test.sample.names))
        ),
    colData=colData
    )

  ds = GenoSet(
    rowRanges=locs,
    assays=list(
        lrr=lrr,
        baf=baf
        ),
    colData=colData
    )

  subset.rows.ds = GenoSet(
    rowRanges=locs[2:3,],
    assays=list(
        lrr=lrr[2:3,],
        baf=baf[2:3,]
        ),
    colData=colData
    )

  subset.cols.ds = GenoSet(
    rowRanges=locs,
    assays=list(
        lrr=matrix(11:30,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3])),
        baf=matrix(41:60,nrow=10,ncol=2,dimnames=list(probe.names,test.sample.names[2:3]))
        ),
    colData=colData[2:3,]
    )

  gene.gr = GRanges(ranges=IRanges(start=2:3,width=1),seqnames=c("chr1","chr1"))

  # Subsetting whole object
  expect_equal( ds[ ,2:3], subset.cols.ds, check.attributes=FALSE)
  expect_equal( ds[ 2:3, ], subset.rows.ds, check.attributes=FALSE)
  expect_equal( ds[ gene.gr, ], subset.rows.ds, check.attributes=FALSE)

  # Subsetting assayData / extracting
  expect_equal( ds[ 5, 3, "baf"], assay(ds,"baf")[5,3])
  expect_equal( ds[ , , "lrr"], assay(ds,"lrr"), label = "Extract whole matrix" )

  # Test subsetting by location
  expect_equal( test.ds[test.gr,], expected.ds, check.attributes=FALSE)
  expect_equal( test.ds[8:10,], expected.ds, check.attributes=FALSE)
  expect_equal( test.ds[ chrIndices(test.ds,"chr3"), ], chr3.ds, check.attributes=FALSE)

  # Replace
  ds = GenoSet(
    rowRanges=locs,
    assays=list(
        lrr=lrr,
        baf=baf
        ),
    colData=colData
    )

  ds[,,"baf"] = ds[,,"lrr"]
  expect_equal(ds[,,"baf"],ds[,,"lrr"], label = "Replace whole element")
  bad.names.lrr = ds[,,"lrr"]
  rownames(bad.names.lrr)[1] = "FOO"
  colnames(bad.names.lrr)[1] = "FOO"
  expect_error({ds[,,"baf"] = bad.names.lrr}, label = "Incoming ad element must have dimnames that matches genoset.",silent=TRUE)
  lrr.mat = ds[,,"lrr"]
  lrr.mat[1:2,1:2] = 5
  ds[1:2,1:2,"lrr"] = 5
  expect_equal(lrr.mat,ds[,,"lrr"], label = "Replace partial matrix with integer indices")
  lrr.mat[6:8,2] = 3
  ds[rowRanges(ds)[6:8,],2,"lrr"] = 3
  expect_equal(lrr.mat,ds[,,"lrr"], label = "Replace partial matrix with RangedData subsetting of rows")
  ds[,3,"lrr"] = 3
  lrr.mat[,3] = 3
  expect_equal(lrr.mat,ds[,,"lrr"], label = "Replace column")
})

test_that("We can sort the genome", {
  chr.names = c(rep("chr1",3),rep("chr2",3),rep("chr10",4))

  ok.locs = GRanges( ranges = IRanges(start=1:10,width=1,names=paste("p",1:10,sep="")), seqnames=factor(chr.names,levels=c("chr1","chr2","chr10")))
  expect_true( isGenomeOrder(ok.locs), label = "Good locs" )

  bad.locs = GRanges( ranges = IRanges(start=c(2,3,1,4,6,5,10:7),width=1,names=paste("p",c(2,3,1,4,6,5,10:7),sep="")), seqnames=factor(chr.names,levels=c("chr1","chr2","chr10")))
  bad.locs.bad.chr = GRanges( ranges = IRanges(start=c(2,3,1,4,6,5,10:7),width=1,names=paste("p",c(2,3,1,4,6,5,10:7),sep="")), seqnames=factor(chr.names,levels=c("chr2","chr1","chr10")))
  expect_true( ! isGenomeOrder(bad.locs, strict=TRUE), label = "Bad within chr, OK chr levels, fail")

  bad.locs2 = GRanges( ranges = IRanges(start=c(1:10),width=1,names=paste("p",c(2,3,1,4,6,5,10:7),sep="")), seqnames=factor(chr.names,levels=c("chr2","chr1","chr10")))
  expect_true( isGenomeOrder(bad.locs2,strict=FALSE) )
  expect_true( ! isGenomeOrder(bad.locs2,strict=TRUE) )

  good.ds = GenoSet(
    rowRanges=ok.locs,
    assays=list(
        cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(rownames(ok.locs),test.sample.names))
        ),
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
    )
  bad.ds = good.ds[ c(6,5,4,3,2,1,10:7),]
  bad.ds.bad.chrs = GenoSet(
    rowRanges=bad.locs.bad.chr,
    assays=list(
        cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(rownames(ok.locs),test.sample.names))[c(4,6,5,2,3,1,10:7),]
        ),
    colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
    )
  expect_true(isGenomeOrder(good.ds))
  expect_true(!isGenomeOrder(bad.ds))
  expect_equal( good.ds, toGenomeOrder(bad.ds,strict=TRUE), check.attributes=FALSE, label = "GenoSet disordered within chrs" )
  expect_equal( good.ds, toGenomeOrder(bad.ds.bad.chrs,strict=TRUE), check.attributes=FALSE, label = "GenoSet disordered within chrs, disordered chrs" )

  gr1 = GRanges(ranges=IRanges(start=c(9,1,5,4,6,2),width=1,names=LETTERS[c(9,1,5,4,6,2)]),seqnames=Rle(factor(c("A","B","C","C","B","A"),levels=c("A","C","B"))))
  gr2 = GRanges(ranges=IRanges(start=c(2,9,4,5,1,6),width=1,names=LETTERS[c(2,9,4,5,1,6)]),seqnames=Rle(factor(c("A","A","C","C","B","B"),levels=c("A","C","B"))))
  gr3 = GRanges(ranges=IRanges(start=c(2,9,1,6,4,5),width=1,names=LETTERS[c(2,9,1,6,4,5)]),seqnames=Rle(factor(c("A","A","B","B","C","C"),levels=c("A","B","C"))))
  expect_identical(toGenomeOrder(gr1,strict=FALSE),gr2, label = "GRanges with mis-ordered chromosomes, without strict")
  expect_identical(toGenomeOrder(gr1,strict=TRUE),gr3, label = "GRanges with mis-ordered chromosomes, with strict")
  expect_true(isGenomeOrder(gr2,strict=FALSE))
  expect_true(isGenomeOrder(gr3,strict=TRUE))
  expect_true(!isGenomeOrder(gr2,strict=TRUE))
  expect_true(!isGenomeOrder(gr1,strict=TRUE), label = "Not in blocks by chromosome, strict")
  expect_true(!isGenomeOrder(gr1,strict=FALSE), label = "Not in blocks by chromosome, strict")
})
