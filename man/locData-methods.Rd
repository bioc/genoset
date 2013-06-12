\docType{methods}
\name{locData}
\alias{locData}
\alias{locData-methods}
\alias{locData,GenoSet-method}
\alias{locData<-}
\alias{locData<--methods}
\alias{locData<-,GenoSet,RangedDataOrGenomicRanges-method}
\title{Access the feature genome position info}
\arguments{
  \item{object}{GenoSet}

  \item{value}{RangedData describing features}
}
\value{
  A GenoSet object
}
\description{
  The position information for each probe/feature is stored
  as an IRanges RangedData object. The locData functions
  allow this data to be accessed or re-set.
}
\section{Methods}{
  \describe{ \item{\code{signature(object = "GenoSet")}}{
  Get location data. } \item{\code{signature(object =
  "GenoSet", value = "RangedData")}}{ Set location data. }}
}
\examples{
data(genoset)
rd = locData(genoset.ds)
locData(genoset.ds) = rd
}
\author{
  Peter M. Haverty
}

