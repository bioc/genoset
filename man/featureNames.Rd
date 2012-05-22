\name{featureNames}
\alias{featureNames}
\alias{featureNames,GRanges-method}
\alias{featureNames,RangedData-method}
\title{Get rownames from RangedData, GRanges, or eSet}
\arguments{
  \item{object}{GRanges, RangedData, or GenoSet}
}
\value{
  character vector with names rows/features
}
\description{
  Get rownames from RangedData, GRanges, or eSet
}
\examples{
data(genoset)
  head(featureNames(locData.rd))
  head(featureNames(as(locData.rd,"GRanges")))
  head(featureNames(cn.ds))
}
\author{
  Peter M. Haverty
}

