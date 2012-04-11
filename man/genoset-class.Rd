\name{genoset-class}
\docType{class}
\alias{class:GenoSet}
\alias{genoset-class}
\alias{GenoSet-class}

\section{Slots}{\describe{\item{\code{locData}:}{(\code{\link{RangedData}}) Contains a RangedData that holds probe locations}}}
\section{Extends}{\code{\linkS4class{eSet}}}

% GenoSet methods
\alias{show,GenoSet-method}
\alias{featureNames<-,GenoSet-method}

% RangedData API added to GenoSet
\alias{universe,GenoSet-method}
\alias{universe<-,GenoSet-method}
\alias{space,GenoSet-method}
\alias{start,GenoSet-method}
\alias{end,GenoSet-method}
\alias{names,GenoSet-method}
\alias{ranges,GenoSet-method}
\alias{elementLengths,GenoSet-method}
\alias{elementLengths,GRanges-method}

% Subsetting
\alias{[,GenoSet,ANY,ANY,ANY-method}
\alias{[,GenoSet,RangedData,ANY,ANY-method}
\alias{[,GenoSet,RangesList,ANY,ANY-method}
\alias{[,GenoSet,character,ANY,ANY-method}

% Genome info from GenoSet or RangedData
\alias{chr,GenoSet-method}
\alias{chr,RangedData-method}
\alias{chr,GRanges-method}
\alias{pos,RangedDataOrGenoSetOrGRanges-method}
\alias{uniqueChrs,RangedDataOrGenoSet-method}
\alias{orderedChrs,RangedDataOrGenoSet-method}
\alias{chrIndices,RangedDataOrGenoSetOrGRanges-method}
\alias{chrIndices,RangedDataOrGenoSetOrGRanges-method}
\alias{chrInfo,RangedDataOrGenoSetOrGRanges-method}
\alias{genoPos,RangedDataOrGenoSet-method}
\alias{chrNames,GenoSet-method}
\alias{chrNames,RangedData-method}
\alias{chrNames,GRanges-method}

\title{GenoSet class}

\description{

The genoset package offers an extension of the BioConductor eSet object for genome
arrays. The package offers three classes. The first class is the
\code{GenoSet} class which can hold an arbitrary number of equal-sized
matrices in its assayData slot. The principal addition of the GenoSet
class is a \code{locData} slot that holds a RangedData object from the
IRanges package. The locData slot allows for quick subsetting by genome position.

Two classes extend GenoSet: CNSet and BAFSet.  CNSet is the basic copy number object. It
keeps its data in the \code{cn} slot, similar to the \code{exprs} slot of the ExpressionSet.
BAFSet is intended to store \code{LRR} or Log-R Ratio and \code{BAF} or B-Allele Frequency data for
SNP arrays.  LRR and BAF come from the terms coined by Illumina. LRR is copynumber data
processed on a per-snp basis to remove some variability using the expected log-ratio of normal
samples with the same genotype. BAF represents the fraction of signal coming from the ``B''
allele, relative to the ``A'' allele, where A and B are arbitrarily assigned. BAF has the
expected value of 0 or 1 for HOM alleles and 0.5 for HET alelles.  Deviation from these
expected values can be interpreted as Allelic Imbalance, which is a sign of gain, loss, or
copy-neutral LOH.
  }


\author{Peter M. Haverty}

\seealso{
  \link{bafset-class},
  \link{cnset-class}
}
  
\examples{
## Creating a GenoSet
test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]
gs = GenoSet(
   locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
   cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
   pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
   annotation="SNP6"
)
}

\keyword{methods}
\keyword{classes}
