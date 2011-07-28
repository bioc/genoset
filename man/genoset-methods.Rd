\name{genoset-methods}
\alias{space,-method}
\alias{space}
\title{Get space factor for GenoSet...}
\description{Get space factor for GenoSet}
\details{\code{space,-method}: locData slot holds a RangedData, which keeps the chromsome of each
feature in a factor names 'space'.

\code{start,-method}: locData slot holds a RangedData.


\code{end,-method}: locData slot holds a RangedData.


\code{names,-method}: Get chromosome names, which are the names of the locData slot.


\code{ranges,-method}: Get ranges from locData slot


\code{elementLengths,-method}: Get elementLengths from locData slot


\code{loadGC-methods}: Local GC content  can be used to remove GC artifacts from copynumber data
see Diskin, 2008). GC% column will be added to the feature data.  The dataset
may be truncated to remove positions without GC information.  GC data are
accessible with locData().


}
\value{\code{space,-method}: factor

\code{start,-method}: integer


\code{end,-method}: integer


\code{names,-method}: character


\code{ranges,-method}: character


\code{elementLengths,-method}: character


\code{loadGC-methods}: An updated object, with GC percentage information added to the locData slot.


}
\author{Peter M. Haverty}
\arguments{\item{x}{GenoSet}
\item{i}{character, RangedData, RangesList, logical, integer}
\item{j}{character, RangedData, RangesList, logical, integer}
\item{drop}{logical drop levels of space factor?}
\item{...}{additional subsetting args}
\item{ds}{A GenoSet object or derivative}
\item{expand}{numeric, expand each feature location by this many bases on each side}
\item{bsgenome,}{sequence db object from BSgenome (e.g. Hsapiens)x}
}
\examples{data(genoset)
space(genoset.ds)
start(genoset.ds)
end(genoset.ds)
names(genoset.ds)
ranges(genoset.ds) # Returns a RangesList
elementLengths(genoset.ds) # Returns the number of probes per chromosomedata(genoset)
genoset.ds[1:5,2:3]  # first five probes and samples 2 and 3
genoset.ds[ , "K"]  # Sample called K
rd = RangedData(ranges=IRanges(start=seq(from=15e6,by=1e6,length=7),width=1),names=letters[8:14],space=rep("chr17",7))
genoset.ds[ rd, "K" ]  # sample K and probes overlapping those in rd, which overlap specifed ranges on chr17
genoset.ds[[ "chr8" ]]    # All samples and probes for chromosome 8}
\alias{start,-method}
\alias{start}
\alias{end,-method}
\alias{end}
\alias{names,-method}
\alias{names}
\alias{ranges,-method}
\alias{ranges}
\alias{elementLengths,-method}
\alias{elementLengths}
\alias{[,GenoSet,ANY,ANY-method}
\alias{[}
\usage{\S4method{[}{GenoSet,ANY,ANY}(x, i, j, ..., drop=FALSE)
\S4method{[}{GenoSet,character,ANY}(x, i, j, ..., drop=FALSE)
\S4method{[}{GenoSet,RangedData,ANY}(x, i, j, ..., drop=FALSE)
\S4method{[}{GenoSet,RangesList,ANY}(x, i, j, ..., drop=FALSE)
\S4method{[[}{GenoSet,character}(x, i, ..., drop=FALSE)
}
\alias{[,GenoSet,character,ANY-method}
\alias{[}
\alias{[,GenoSet,RangedData,ANY-method}
\alias{[}
\alias{[,GenoSet,RangesList,ANY-method}
\alias{[}
\alias{[[,GenoSet,character-method}
\alias{[[}
\alias{loadGC-methods}
\alias{loadGC}
\alias{loadGC,-method}
\alias{loadGC}
\alias{loadGC,-method}
\alias{loadGC}

