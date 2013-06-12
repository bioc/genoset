\name{featureNames}
\alias{featureNames}
\alias{featureNames,GenoSet-method}
\alias{featureNames,GRanges-method}
\alias{featureNames,RangedData-method}
\title{Get rownames from RangedData, GRanges, or GenoSet}
\arguments{
  \item{object}{GRanges, RangedData, or GenoSet}
}
\value{
  character vector with names rows/features
}
\description{
  Get rownames from RangedData, GRanges, or GenoSet
}
\examples{
data(genoset)
  head(featureNames(locData.gr))
  head(featureNames(genoset.ds))
}
\author{
  Peter M. Haverty
}

