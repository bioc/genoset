\docType{methods}
\name{featureNames,GenoSet-method}
\alias{featureNames,GenoSet-method}
\title{Get rownames from RangedData, GRanges, or GenoSet}
\usage{
\S4method{featureNames}{GenoSet}(object)
}
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
\author{
Peter M. Haverty
}

