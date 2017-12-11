library(RUnit)
library(genoset)

sample.names = LETTERS[11:13]
probe.names = letters[1:10]
locs = GRanges(ranges=IRanges(start=c(1,3,5,7,4,6,2,4,6,8),width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))

checkEquals(
    chrIndices(locs),
    matrix( c(1,5,7,4,6,10,0,4,6), ncol=3, dimnames=list(c("chr1","chr3","chrX"),c("first","last","offset")) )
    )
           
                                                          
