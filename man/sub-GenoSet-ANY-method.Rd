\docType{methods}
\name{[,GenoSet,ANY-method}
\alias{[,GenoSet,ANY,ANY,ANY-method}
\alias{[,GenoSet,ANY-method}
\title{Subset a GenoSet}
\usage{
\S4method{[}{GenoSet,ANY}(x, i, j, k, ..., drop = FALSE)
}
\arguments{
  \item{x}{GenoSet}

  \item{i}{character, RangedData, logical, integer}

  \item{j}{character, RangedData, logical, integer}

  \item{k}{character or integer}

  \item{drop}{logical drop levels of space factor?}

  \item{...}{additional subsetting args}
}
\description{
Subset a GenoSet
}
\examples{
data(genoset)
  genoset.ds[1:5,2:3]  # first five probes and samples 2 and 3
  genoset.ds[ , "K"]  # Sample called K
  rd = RangedData(ranges=IRanges(start=seq(from=15e6,by=1e6,length=7),width=1),names=letters[8:14],space=rep("chr17",7))
  genoset.ds[ rd, "K" ]  # sample K and probes overlapping those in rd, which overlap specifed ranges on chr17
}

