\name{bafset-class}
\docType{class}
\alias{class:BAFSet}
\alias{bafset-class}
\alias{BAFSet-class}

\section{Extends}{\code{\linkS4class{GenoSet}}}

\alias{baf,BAFSet-method}
\alias{lrr,BAFSet-method}
\alias{baf<-,BAFSet,matrix-method}
\alias{lrr<-,BAFSet,matrix-method}
  
% Coercion:
\alias{coerce,CNSet,BAFSet-method}
\alias{coerce,ExpressionSet,BAFSet-method}

\title{BAFSet class}
\description{A BAFSet is and extension of GenoSet that requires 'baf'
  and 'lrr' assayData element
}

\author{Peter M. Haverty}

\seealso{
  \link{bafset-class},
  \link{cnset-class}
}

\examples{
## Creating a BAFSet
test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]
locData.rd = RangedData(ranges=IRanges(start=c(1,4,3,2,5:10),width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18")
bs = BAFSet(
  locData=locData.rd,
  lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
  baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
  pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
  annotation="SNP6"
)
}

\keyword{methods}
\keyword{classes}
