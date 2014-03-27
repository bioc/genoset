\docType{methods}
\name{featureNames,GenoSet-method}
\alias{featureNames,GRanges-method}
\alias{featureNames,GenoSet-method}
\alias{featureNames,RangedData-method}
\alias{featureNames<-,GRanges-method}
\alias{featureNames<-,GenoSet-method}
\alias{featureNames<-,RangedData-method}
\alias{rownames,GRanges-method}
\alias{rownames,GenoSet-method}
\alias{rownames<-,GRanges-method}
\title{Get rownames from RangedData, GRanges, or GenoSet}
\usage{
\S4method{featureNames}{GenoSet}(object)

\S4method{featureNames}{GRanges}(object)

\S4method{featureNames}{RangedData}(object)

\S4method{rownames}{GRanges}(x)

\S4method{rownames}{GenoSet}(x)

\S4method{featureNames}{GenoSet}(object) <- value

\S4method{featureNames}{GRanges}(object) <- value

\S4method{featureNames}{RangedData}(object) <- value

\S4method{rownames}{GRanges}(x) <- value
}
\arguments{
  \item{object}{GRanges, RangedData, or GenoSet}

  \item{value}{Character vector containing feature or
  sample names.}

  \item{x}{ A matrix-like object.  }
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

