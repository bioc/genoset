\name{chr}
\alias{chr}
\alias{chr,GenoSet-method}
\alias{chr,GRanges-method}
\alias{chr,RangedData-method}
\title{Look up chromosome for each feature}
\arguments{
  \item{object}{GRanges, RangedData or GenoSet}
}
\value{
  character vector of chromosome positions for each feature
}
\description{
  Chromosome name for each feature
}
\details{
  Get chromosome name for each feature.  Returns character,
  not the factor 'space'.
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
  chr(gs)  # c("chr1","chr1","chr1","chr1","chr3","chr3","chrX","chrX","chrX","chrX")
  chr(locData(gs))  # The same
}
\author{
  Peter Haverty
}

