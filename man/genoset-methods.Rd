\name{universe}
\alias{[}
\alias{[,GenoSet,ANY-method}
\alias{[,GenoSet,ANY,ANY,ANY-method}
\alias{[,GenoSet,character-method}
\alias{[,GenoSet,character,ANY,ANY-method}
\alias{[,GenoSet,RangedDataOrRangesListOrGRanges-method}
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

  integer

  integer

  integer

  character

  character

  factor

  character
}
\description{
  Genome universe for locData

  Set genome universe

  Get start of location for each feature

  Get end of location for each feature

  Get width of location for each feature

  Get chromosome names

  Get ranges from locData slot

  locData slot holds a RangedData, which keeps the
  chromosome of each feature in a factor names 'space'. The
  ranges method on a GenoSet is depricated. Please use
  space(locData(x)) or seqnames(locData(x)) as appropriate
  for RangedData or GRanges.

  Get elementLengths from locData slot
}
\details{
  The genome positions of the features in locData. The UCSC
  notation (e.g. hg18, hg19, etc.) should be used.

  Get chromosome names, which are the names of the locData
  slot. The names method on a GenoSet is depricated. Please
  use chrNames.

  Get ranges from locData slot. The ranges method on a
  GenoSet is depricated. Please use ranges(locData(x)).

  Get elementLengths from locData slot
}
\examples{
data(genoset)
  universe(locData.rd)
  universe(locData.rd) = "hg19"
data(genoset)
chr(genoset.ds)
start(genoset.ds)
end(genoset.ds)
chrNames(genoset.ds)
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

  Peter Haverty

  Peter Haverty

  Peter M. Haverty

  Peter Haverty
}

