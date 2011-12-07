test.assayData <-function() {
  
  test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  
  ds = BAFSet(
    locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
    )
  bm.ds = convertToBigMatrix(ds,path=tempdir())

  assayDataElement(bm.ds,"baf") <- ds[,,"baf"]
  checkEquals( ds[,,"baf"], bm.ds[,,"baf"], "Replace writeable big.matrix OK" )
  Sys.chmod(attr(bm.ds[,,"lrr"],"desc"),"0444")
  checkException( assayDataElement(bm.ds,"lrr") <- ds[,,"lrr"], silent=TRUE, "Replacing non-writeable big.matrix not OK" )
  checkEquals( bm.ds[,,"lrr"][,], ds[,,"lrr"], "Failing to modify lrr in above test should not modify lrr" )
  rdata.file = file.path(tempdir(),"bm.RData")
  save(bm.ds,file=rdata.file)
  loaded.ds = readGenoSet(rdata.file)
  checkTrue( is.big.matrix(assayDataElement(loaded.ds,"lrr")) & !is.nil(assayDataElement(loaded.ds,"lrr")@address))
}
