\name{BAFSet-class}
\Rdversion{1.1}
\docType{class}
\alias{BAFSet-class}

\title{Class \code{"BAFSet"}}
\description{A BAFSet is and extension of GenoSet that requires 'baf' and 'lrr' assayData element}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("BAFSet", assayData, phenoData, featureData, experimentData, annotation, protocolData, ...)}.
However, as per BioConductor standard practice the object creation function \code{BAFSet} is recommended.
}
\section{Slots}{
  \describe{
    \item{\code{locData}:}{Object of class \code{"RangedData"} Feature
      locations on the genome }
    \item{\code{assayData}:}{Object of class \code{"AssayData"} ~~ }
    \item{\code{phenoData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{featureData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{experimentData}:}{Object of class \code{"MIAxE"} ~~ }
    \item{\code{annotation}:}{Object of class \code{"character"} ~~ }
    \item{\code{protocolData}:}{Object of class \code{"AnnotatedDataFrame"} ~~ }
    \item{\code{.__classVersion__}:}{Object of class \code{"Versions"} ~~ }
  }
}
\section{Extends}{
Class \code{"\linkS4class{GenoSet}"}, directly.
Class \code{"\linkS4class{eSet}"}, by class "GenoSet", distance 2.
Class \code{"\linkS4class{RangedDataOrGenoSet}"}, by class "GenoSet", distance 2.
Class \code{"\linkS4class{RangedDataOrGenoSetOrGRanges}"}, by class "GenoSet", distance 2.
Class \code{"\linkS4class{VersionedBiobase}"}, by class "GenoSet", distance 3.
Class \code{"\linkS4class{Versioned}"}, by class "GenoSet", distance 4.
}
\section{Methods}{
  \describe{
    \item{baf}{\code{signature(object = "BAFSet")}: Getter for 'baf' assayDataElement }
    \item{baf<-}{\code{signature(object = "BAFSet", value = "matrix")}: setter for 'baf' assayDataElement }
    \item{genoPlot}{\code{signature(x = "BAFSet", y = "ANY")}: Plot data
      along the genome. Defaults to 'lrr' assayDataElement. }x
    \item{lrr}{\code{signature(object = "BAFSet")}: Getter for 'lrr' assayDataElement }
    \item{lrr<-}{\code{signature(object = "BAFSet", value = "matrix")}: Setter for 'lrr' assayDataElement }
	 }
}
\author{
Peter M. Haverty \email{phaverty@gene.com}
}
\seealso{
  \code{\link{BAFSet}}, \code{\linkS4class{CNSet}}, \code{\linkS4class{GenoSet}}
}
\examples{
showClass("BAFSet")
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
\keyword{classes}
