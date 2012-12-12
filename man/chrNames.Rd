\name{chrNames}
\alias{chrNames}
\alias{chrNames,GenoSet-method}
\alias{chrNames,GRanges-method}
\alias{chrNames,RangedData-method}
\alias{chrNames<-}
\alias{chrNames<-,GenoSet-method}
\alias{chrNames<-,GRanges-method}
\alias{chrNames<-,RangedData-method}
\title{Get list of unique chromosome names}
\arguments{
  \item{object}{RangedData or GenoSet}
}
\value{
  character vector with names of chromosomes
}
\description{
  Get list of unique chromosome names
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
  chrNames(gs) # c("chr1","chr3","chrX")
  chrNames(locData(gs))  # The same
  chrNames(gs) = sub("^chr","",chrNames(gs))
}
\author{
  Peter M. Haverty
}

