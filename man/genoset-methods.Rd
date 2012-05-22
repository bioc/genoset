\name{universe}
\alias{[}
\alias{[,GenoSet,ANY,ANY,ANY-method}
\alias{[,GenoSet,character,ANY,ANY-method}
\alias{[,GenoSet,RangedDataOrRangesListOrGRanges,ANY,ANY-method}
\alias{[<-}
\alias{[<-,GenoSet,ANY,ANY,ANY-method}
\alias{elementLengths}
\alias{elementLengths,GenoSet-method}
\alias{elementLengths,GRanges-method}
\alias{end}
\alias{end,GenoSet-method}
\alias{names}
\alias{names,GenoSet-method}
\alias{ranges}
\alias{ranges,GenoSet-method}
\alias{show}
\alias{show,GenoSet-method}
\alias{space}
\alias{space,GenoSet-method}
\alias{start}
\alias{start,GenoSet-method}
\alias{universe}
\alias{universe,GenoSet-method}
\alias{universe<-}
\alias{universe<-,GenoSet-method}
\alias{width}
\title{Get and set the genome universe annotation.}
\arguments{
  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{value}{character, new universe string, e.g. hg19}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{i}{character, RangedData, RangesList, logical,
  integer}

  \item{j}{character, RangedData, RangesList, logical,
  integer}

  \item{k}{character or integer}

  \item{drop}{logical drop levels of space factor?}

  \item{...}{additional subsetting args}
}
\value{
  character, e.g. hg19

  A GenoSet object

  factor

  integer

  integer

  integer

  character

  character

  character
}
\description{
  Genome universe for locData

  Set genome universe

  locData slot holds a RangedData, which keeps the
  chromosome of each feature in a factor names 'space'.

  locData slot holds a RangedData.

  locData slot holds a RangedData.

  locData slot holds a RangedData.

  Get chromosome names

  Get ranges from locData slot

  Get elementLengths from locData slot
}
\details{
  The genome positions of the features in locData. The UCSC
  notation (e.g. hg18, hg19, etc.) should be used.

  Get chromosome names, which are the names of the locData
  slot.

  Get ranges from locData slot

  Get elementLengths from locData slot
}
\examples{
data(genoset)
  universe(genoset.ds)
  universe(genoset.ds) = "hg19"
data(genoset)
space(genoset.ds)
start(genoset.ds)
end(genoset.ds)
chrNames(genoset.ds)
ranges(genoset.ds) # Returns a RangesList
elementLengths(genoset.ds) # Returns the number of probes per chromosome
data(genoset)
  genoset.ds[1:5,2:3]  # first five probes and samples 2 and 3
  genoset.ds[ , "K"]  # Sample called K
  rd = RangedData(ranges=IRanges(start=seq(from=15e6,by=1e6,length=7),width=1),names=letters[8:14],space=rep("chr17",7))
  genoset.ds[ rd, "K" ]  # sample K and probes overlapping those in rd, which overlap specifed ranges on chr17
}
\author{
  Peter M. Haverty

  Peter Haverty

  Peter M. Haverty

  Peter M. Haverty

  Peter M. Haverty

  Peter M. Haverty

  Peter Haverty

  Peter Haverty

  Peter Haverty
}

