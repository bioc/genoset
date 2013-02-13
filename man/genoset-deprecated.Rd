\name{genoset-deprecated}
\alias{genoset-deprecated}
\title{Deprecated genoset features}
\description{
  The CNSet and BAFSet classes have been deprecated.  They
  only really added getter/setter methods for specific
  assayDataElements, so they are now redundant with the
  preferred method of using the assayDataElement name as
  the third argument to bracket, e.g. \code{x[i, j,
  "lrr"]}. Accordingly \code{BAFSet.to.ExpressionSets} is
  also deprecated.
}
\details{
  Additionally, names, ranges, and space on a GenoSet are
  also deprecated. In an effort to make a consistent API
  for either RangedData or GRanges in the locData slot, we
  recommend using \code{chrNames} for \code{names} and
  \code{chr} for \code{space}.
}

