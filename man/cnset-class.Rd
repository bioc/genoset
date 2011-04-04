\name{cnset-class}
\docType{class}
\alias{class:CNSet}
\alias{cnset-class}
\alias{CNSet-class}

\section{Extends}{\code{\linkS4class{GenoSet}}}

\alias{cn,CNSet-method}
\alias{cn<-,CNSet,matrix-method}
\alias{genoPlot,CNSet,ANY-method}

% Coercion:
\alias{coerce,ExpressionSet,CNSet-method}

\title{CNSet class}
\description{A CNSet is an extension of GenoSet that requires a 'cn'
  assayData element.
}

\author{Peter M. Haverty}

\seealso{
  \link{bafset-class},
  \link{cnset-class}
}

\examples{
test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]
cn.ds = CNSet(
   locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
   cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
   pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
   annotation="SNP6"
   )
}

\keyword{methods}
\keyword{classes}
