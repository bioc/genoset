\name{chrInfo}
\alias{chrInfo-methods}
\alias{chrInfo}
\title{Chromosome Information}
\description{Get chromosome start and stop positions}
\details{\code{chrInfo-methods}: Provides a matrix of start, stop and offset, in base numbers for each chromosome.

}
\value{\code{chrInfo-methods}: list with start and stop position, by ordered chr

}
\author{Peter Haverty}
\arguments{\item{object}{A GenoSet object or similar}
}
\examples{data(genoset)
chrInfo(genoset.ds)
chrInfo(locData(genoset.ds))  # The same}
\alias{chrInfo,-method}
\alias{chrInfo}

