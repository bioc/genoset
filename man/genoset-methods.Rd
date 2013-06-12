\name{universe}
\alias{[}
\alias{[,GenoSet,ANY-method}
\alias{[,GenoSet,ANY,ANY,ANY-method}
\alias{[,GenoSet,character-method}
\alias{[,GenoSet,character,ANY,ANY-method}
\alias{[,GenoSet,RangedDataOrGenomicRanges-method}
\alias{[,GenoSet,RangedDataOrGenomicRanges,ANY,ANY-method}
\alias{[<-}
\alias{[<-,GenoSet,ANY,ANY,ANY-method}
\alias{dim}
\alias{dim,GenoSet-method}
\alias{elementLengths}
\alias{elementLengths,GenoSet-method}
\alias{elementLengths,GRanges-method}
\alias{end}
\alias{end,GenoSet-method}
\alias{names}
\alias{names,GenoSet-method}
\alias{nrow}
\alias{nrow,GRanges-method}
\alias{start}
\alias{start,GenoSet-method}
\alias{universe}
\alias{universe,GenoSet-method}
\alias{universe,GRanges-method}
\alias{universe<-}
\alias{universe<-,GenoSet-method}
\alias{universe<-,GRanges-method}
\alias{width}
\title{Get and set the genome universe annotation.}
\arguments{
  \item{x}{GenoSet or GRanges}

  \item{x}{GenoSet or GRanges}

  \item{value}{character, new universe string, e.g. hg19}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{i}{character, RangedData, logical, integer}

  \item{j}{character, RangedData, logical, integer}

  \item{k}{character or integer}

  \item{drop}{logical drop levels of space factor?}

  \item{...}{additional subsetting args}
}
\value{
  character, e.g. hg19

  updated copy of x

  integer

  integer

  integer

  character

  character
}
\description{
  Genome universe for locData

  Set genome universe

  Get start of location for each feature

  Get end of location for each feature

  Get width of location for each feature

  Get names of data matrices. For the time being, this is
  \code{assayDataElementNames}. This function used to do
  \code{chrNames}.

  Get elementLengths from locData slot
}
\details{
  The genome positions of the features in locData. The UCSC
  notation (e.g. hg18, hg19, etc.) should be used. For a
  GRanges, the first value is returned if there are
  multiple.

  Get elementLengths from locData slot
}
\examples{
data(genoset)
  universe(locData.gr)
  universe(locData.gr) = "hg19"
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
}

