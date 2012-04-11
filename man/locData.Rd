\name{locData}
\alias{locData}
\alias{locData,GenoSet-method}
\alias{locData<-}
\alias{locData<-,GenoSet,RangedData-method}
\title{Get and set probe set info}
\arguments{
  \item{object}{GenoSet}

  \item{object}{A GenoSet object}

  \item{object}{GenoSet}

  \item{value}{RangedData describing features}
}
\value{
  A GenoSet object
}
\description{
  Access the feature genome position info

  Set locData
}
\details{
  The position information for each probe/feature is stored
  as an IRanges RangedData object. The locData functions
  allow this data to be accessed or re-set.

  Set locData
}
\examples{
data(genoset)
  rd = locData(genoset.ds)
  locData(genoset.ds) = rd
}
\author{
  Peter M. Haverty

  Peter Haverty
}

