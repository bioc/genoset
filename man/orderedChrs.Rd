\name{orderedChrs}
\alias{orderedChrs}
\title{Get chromosome names in genome order}
\arguments{
  \item{object}{GenoSet or RangedData}
}
\value{
  character vector with chrs in genome order
}
\description{
  Get chromosome names from locData data in a GenoSet.
  Order numerically, for numeric chromosomes, then
  lexically for the rest.
}
\examples{
test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  gs = GenoSet(
     locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chrX",2),rep("chr3",4)),universe="hg18"),
     cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
     pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
     annotation="SNP6"
  )
  orderedChrs(gs) # c("chr1","chr3","chrX")
  orderedChrs(locData(gs))  # The same
}
\author{
  Peter M. Haverty
}

