\name{boundingIndicesByChr}
\alias{boundingIndicesByChr}
\title{Find indices of features bounding a set of chromosome ranges/genes, across chromosomes}
\usage{
  boundingIndicesByChr(query, subject)
}
\arguments{
  \item{query}{RangedData or GRanges (with genes). GRanges
  coerced to RangedData for now.}

  \item{subject}{RangedData}
}
\value{
  integer matrix with two columns corresponding to indices
  on left and right bound of queries in subject
}
\description{
  Finds subject ranges corresponding to a set of genes
  (query ranges), taking chromosome into account.
  Specifically, this function will find the indices of the
  features (first and last) bounding the ends of a
  range/gene (start and stop) such that first <= start <
  stop <= last. Equality is necessary so that multiple
  conversions between indices and genomic positions will
  not expand with each conversion. Ranges/genes that are
  outside the range of feature positions will be given the
  indices of the corresponding first or last index on that
  chromosome, rather than 0 or n + 1 so that genes can
  always be connected to some data. Checking the left and
  right bound for equality will tell you when a query is
  off the end of a chromosome.
}
\details{
  This function uses some tricks from findIntervals, where
  is for k queries and n features it is O(k * log(n))
  generally and ~O(k) for sorted queries. Therefore will be
  dramatically faster for sets of query genes that are
  sorted by start position within each chromosome. The
  index of the stop position for each gene is found using
  the left bound from the start of the gene reducing the
  search space for the stop position somewhat. This
  function has important differences from boundingIndices2,
  which uses findInterval: boundingIndices does not check
  for NAs or unsorted data in the subject positions. Also,
  the positions are kept as integer, where boundingIndices2
  (and findInterval) convert them to doubles. These three
  once-per-call differences account for much of the speed
  improvement in boundingIndices. These three differences
  are meant for position info coming from GenoSet objects
  and boundingIndices2 is safer for general use.
  boundingIndices works on integer postions and does not
  check that the positions are ordered. The starts and
  stops need not be sorted, but it will be much faster if
  they are.

  This function differs from boundingIndices in that 1. it
  uses both start and end positions for the subject, and 2.
  query and subject start and end positions are processed
  in blocks corresponding to chromosomes.
}
\author{
  Peter M. Haverty \email{phaverty@gene.com}
}
\seealso{
  Other "range summaries": \code{\link{boundingIndices}},
  \code{\link{boundingIndices2}},
  \code{\link{rangeColMeans}},
  \code{\link{rangeSampleMeans}}
}

