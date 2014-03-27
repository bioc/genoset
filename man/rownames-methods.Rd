\name{featureNames}
\alias{featureNames}
\alias{featureNames<-}
\alias{rownames}
\alias{rownames<-}
\title{Get rownames from RangedData, GRanges, or GenoSet}
\arguments{
  \item{object}{GRanges, RangedData, or GenoSet}
}
\value{
  character vector with names rows/features
}
\description{
  Get rownames from RangedData, GRanges, or GenoSet.
}
\examples{
data(genoset)
  head(rownames(locData.gr))
  head(rownames(genoset.ds))
}

