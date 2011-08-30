\name{initGenoSet}
\alias{initGenoSet}
\title{Create a GenoSet or derivative object...}
\usage{initGenoSet(type, locData, pData, annotation="", universe, ...)
}
\description{Create a GenoSet or derivative object}
\details{This function is the preferred method for creating a new GenoSet object. Users are
generally discouraged from calling "new" directly. The "..." argument is for any number of matrices of matching size that will
become part of the assayData slot of the resulting object. This function passes
control to the "genoSet" object which performs argument checking including
dimname matching among relevant slots and sets everything to genome order. Genome
order can be disrupted by "[" calls and will be checked by methods that
require it.}
\value{A GenoSet object or derivative as specified by "type" arg}
\author{Peter M. Haverty}
\arguments{\item{type}{character, the type of object (e.g. GenoSet, BAFSet, CNSet) to be created}
\item{locData}{A RangedData object specifying feature chromosome
locations. Rownames are required to match featureNames.}
\item{pData}{A data frame with rownames matching all data matrices}
\item{annotation}{character, string to specify chip/platform type}
\item{universe}{character, a string to specify the genome universe
for locData}
\item{...}{More matrix or DataFrame objects to include in assayData}
}
\examples{test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]
gs = GenoSet(
locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
annotation="SNP6"
)}

