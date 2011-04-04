\name{genoPos-methods}
\alias{genoPos-methods}
\alias{genoPos}
\title{Convert chromosome positions to positions from start of genome}
\description{Get base positions of features in genome-scale units}
\details{\code{genoPos-methods}: Get base positions of array features in bases counting from the start of the
genome. Chromosomes are ordered numerically, when possible, then lexically.

}
\value{\code{genoPos-methods}: numeric position of each feature in whole genome units, in original order

}
\author{Peter M. Haverty}
\arguments{\item{object}{A GenoSet object or a RangedData object}
}
\examples{data(genoset)
genoPos(genoset.ds)
genoPos(locData(genoset.ds))  # The same}
\alias{genoPos,-method}
\alias{genoPos}

