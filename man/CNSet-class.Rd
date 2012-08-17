\name{CNSet-class}
\Rdversion{1.1}
\docType{class}
\alias{CNSet-class}

\title{Class \code{"CNSet"}}
\description{A CNSet is an extension of GenoSet that requires a 'cn' assayData element.}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CNSet",
  assayData, phenoData, featureData, experimentData, annotation,
  protocolData, ...)}.
However, as per BioConductor standard practice the object creation function \code{CNSet} is recommended.
}
\section{Slots}{
  \describe{
    \item{\code{locData}:}{Object of class \code{"RangedDataOrGRanges"} Feature locations on the genome. }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} From eSet }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} From eSet }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"}  From eSet }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} From eSet }
    \item{\code{annotation}:}{Object of class \code{"character"} From eSet }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} From eSet }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} From eSet }
  }
}
\section{Extends}{
Class \code{"\linkS4class{GenoSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "GenoSet", distance 2.
}
\section{Methods}{
  \describe{
    \item{cn}{\code{signature(object = "CNSet")}: Getter for cn assayDataElement }
    \item{cn<-}{\code{signature(object = "CNSet", value = "matrix")}: Setter for 'cn' assayDataElement }
    \item{genoPlot}{\code{signature(x = "CNSet", y = "ANY")}: Plot data
  along the genome. Defaults to 'cn' assayDataElement }
	 }
}
\author{
Peter M. Haverty \email{phaverty@gene.com}
}
\seealso{
  \code{\link{CNSet}}, \code{\linkS4class{GenoSet}}, \code{\linkS4class{BAFSet}}
}
\examples{
showClass("CNSet")
test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]
cn.ds = CNSet(
   locData=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4))),
   cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
   pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
   annotation="SNP6"
   )
}
\keyword{classes}
