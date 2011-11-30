\name{BAFSet.to.ExpressionSets}
\alias{BAFSet.to.ExpressionSets}
\title{Make a pair of ExpressionSets from a BAFSet}
\usage{
  BAFSet.to.ExpressionSets(bs)
}
\arguments{
  \item{bs}{A BAFset object}
}
\value{
  A list with one ExpressionSet each for the baf and lrr
  data in the BAFSet object
}
\description{
  Often it is convenient to have a more standard
  "ExpressionSet" rather than a BAFSet.  For example, when
  using infrastructure dependent on the ExpressionSet
  slots, like limma or ExpressionSetOnDisk. This will
  create a list of two ExpressionSets, one each for the baf
  and lrr data. To make a single ExpressionSet, with the
  lrr data in the exprs slot and the baf data as an
  additional member of assayData, use the standard coercion
  eset = as(bafset,"ExpressionSet").
}
\examples{
data(genoset)
  eset.list = BAFSet.to.ExpressionSets(baf.ds)
}
\author{
  Peter M. Haverty
}

