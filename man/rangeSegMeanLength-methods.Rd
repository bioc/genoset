\name{rangeSegMeanLength}
\alias{rangeSegMeanLength}
\title{Get segment widths}
\arguments{
  \item{range.gr}{GRanges, genome regions of interest,
  usually genes}

  \item{segs}{data.frame of segments, like from segTable,
  or a list of these}
}
\value{
  named vector of lengths, one per item in range.gr, or a
  range x length(segs) of these if segs is also list-like.
}
\description{
  The width of a genomic segment helps inform us about the
  importance of a copy number value. Focal amplifications
  are more interesting than broad gains, for example. Given
  a range of interesting regions (i.e. genes) this function
  determines all genomics segments covered by each gene and
  returns the average length of the segments covered by
  each gene in each sample. Often only a single segment
  covers a given gene in a given sample.
}
\seealso{
  Other "segmented data": \code{\link{bounds2Rle}},
  \code{\link{runCBS}}, \code{\link{segPairTable}},
  \code{\link{segPairTable}}, \code{\link{segPairTable}},
  \code{\link{segs2Granges}},
  \code{\link{segs2RangedData}}, \code{\link{segs2Rle}},
  \code{\link{segs2RleDataFrame}}, \code{\link{segTable}},
  \code{\link{segTable}}, \code{\link{segTable}}
}

