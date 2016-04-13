\name{GenoSet-class}
\Rdversion{1.1}
\docType{class}
\alias{GenoSet-class}
\alias{GenoSetOrGenomicRanges-class}


\title{Class \code{"GenoSet"}}
\description{
\code{GenoSet} extends \code{RangedSummarizedExperiment} by adding some
additional methods to the API. Examples include subsetting rows with a
GenomicRanges and combining this with access to assays like \code{genoset[i,j,assay]}.
}
\section{Extends}{
Class \code{\linkS4class{RangedSummarizedExperiment}}, directly.
}
\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "GenoSet", i = "ANY", j = "ANY", drop = "ANY")}: ... }
    \item{[}{\code{signature(x = "GenoSet", i = "character", j = "ANY", drop = "ANY")}: ... }
    \item{[<-}{\code{signature(x = "GenoSet", i = "ANY", j = "ANY", value = "ANY")}: ... }
    \item{chr}{\code{signature(object = "GenoSet")}: ... }
    \item{chrNames}{\code{signature(object = "GenoSet")}: ... }
    \item{dim}{\code{signature(object = "GenoSet")}: ... }
    \item{genoPlot}{\code{signature(x = "GenoSet", y = "ANY")}: ... }
    \item{rowRanges}{\code{signature(object = "GenoSet")}: ... }
    \item{names}{\code{signature(x = "GenoSet")}: ... }
    \item{ranges}{\code{signature(x = "GenoSet")}: ... }
    \item{chrInfo}{\code{signature(x = "GenoSet")}: ... }
    \item{chrIndices}{\code{signature(x = "GenoSet")}: ... }
    \item{show}{\code{signature(object = "GenoSet")}: ... }
    \item{toGenomeOrder}{\code{signature(ds = "GenoSet")}: ... }
    \item{isGenomeOrder}{\code{signature(ds = "GenoSet")}: ... }
    \item{assays}{\code{signature(x = "GenoSet")}: ... }
    \item{assay}{\code{signature(x = "GenoSet", i="ANY")}: ... }
    \item{assay<-}{\code{signature(x = "GenoSet", i="ANY",value="ANY")}: ... }
    \item{assayNames}{\code{signature(x = "GenoSet")}: ... }
    \item{colData}{\code{signature(x = "GenoSet")}: ... }
    \item{locData}{\code{signature(x = "GenoSet")}: ... }
    \item{locData<-}{\code{signature(x = "GenoSet",value="GenomicRanges")}: ... }
  }
}

\seealso{
\code{\link{GenoSet}}
}
\examples{
showClass("GenoSet")
test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]
assays=list(matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)))
rowRanges=GRanges(ranges=IRanges(start=1:10,width=1,names=probe.names),seqnames=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)))
colData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5])))
rse=SummarizedExperiment(rowRanges=rowRanges,assays=assays,colData=colData,metadata=metadata)
gs = GenoSet(rowRanges, assays, colData)
}
\keyword{classes}
