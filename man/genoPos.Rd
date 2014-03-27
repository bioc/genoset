\name{genoPos}
\alias{genoPos}
\title{Convert chromosome positions to positions from start of genome}
\usage{
genoPos(object)
}
\arguments{
  \item{object}{A GenoSet object or a RangedData object}
}
\value{
numeric position of each feature in whole genome units, in
original order
}
\description{
Get base positions of features in genome-scale units
}
\details{
Get base positions of array features in bases counting from
the start of the genome. Chromosomes are ordered
numerically, when possible, then lexically.
}
\examples{
data(genoset)
  head(genoPos(genoset.ds))
  head(genoPos(locData(genoset.ds)))  # The same
}

