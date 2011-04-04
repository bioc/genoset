\name{locData}
\alias{locData-methods}
\alias{locData}
\title{Get and set probe set info}
\description{Access the feature genome position info}
\details{\code{locData-methods}: The position information for each probe/feature is stored as an IRanges RangedData object.
The locData functions allow this data to be accessed or re-set.

\code{locData<-,-method}: Set locData


}
\author{Peter M. Haverty}
\arguments{\item{object}{GenoSet}
\item{object}{A GenoSet object}
\item{value}{RangedData describing features}
}
\examples{data(genoset)
rd = locData(genoset.ds)
locData(genoset.ds) = rd}
\alias{locData,-method}
\alias{locData}
\alias{locData<-,-method}
\alias{locData<-}
\value{\code{locData<-,-method}: A GenoSet object

}

