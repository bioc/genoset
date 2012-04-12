\name{pos}
\alias{pos}
\alias{pos,RangedDataOrGenoSetOrGRanges-method}
\title{Positions for features}
\arguments{
  \item{object}{GRanges, RangedData or GenoSet}
}
\value{
  numeric vector of feature positions within a chromosome
}
\description{
  Chromosome position of features
}
\details{
  Get chromosome position of features/ranges. Defined as
  floor of mean of start and end.
}
\examples{
test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  gs = GenoSet(
     locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
     cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
     pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
     annotation="SNP6"
  )
  pos(gs)  # 1:10
  pos(locData(gs))  # The same
}
\author{
  Peter Haverty
}

