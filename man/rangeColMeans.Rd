\name{rangeColMeans}
\alias{rangeColMeans}
\title{Calculate column means for multiple ranges}
\usage{
  rangeColMeans(bounds, x)
}
\arguments{
  \item{bounds}{A two column integer matrix of row indices}

  \item{x}{A numeric matrix with rows corresponding to
  indices in bounds.}
}
\value{
  A numeric matrix or vector, matching the form of x. One
  row for each row in bounds, one col for each col of x and
  appropriate dimnames. If x is a vector, just a vector
  with names from the rownames of bounds.
}
\description{
  Essentially colMeans with a loop, all in a .Call.
  Designed to take a 2-column matrix of row indices,
  bounds, for a matrix, x, and calculate mean for each
  range in each column (or along a single vector). bounds
  matrix need not cover all rows.
}
\author{
  Peter M. Haverty \email{phaverty@gene.com}
}
\seealso{
  Other "range summaries": \code{\link{boundingIndices}},
  \code{\link{boundingIndices2}},
  \code{\link{rangeSampleMeans}}
}

