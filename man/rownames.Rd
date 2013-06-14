\name{featureNames}
\alias{featureNames}
\alias{featureNames,GenoSet-method}
\alias{featureNames,GRanges-method}
\alias{featureNames,RangedData-method}
\alias{featureNames<-}
\alias{featureNames<-,GenoSet-method}
\alias{featureNames<-,GRanges-method}
\alias{featureNames<-,RangedData-method}
\alias{rownames}
\alias{rownames,GenoSet-method}
\alias{rownames,GRanges-method}
\alias{rownames<-}
\alias{rownames<-,GRanges-method}
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
\author{
  Peter M. Haverty
}

