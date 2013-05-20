\name{segs2RleDataFrame}
\alias{segs2RleDataFrame}
\title{CBS segments to probe matrix}
\usage{
  segs2RleDataFrame(seg.list, locs)
}
\arguments{
  \item{seg.list}{list, list of data frames, one per
  sample, each is result from CBS}

  \item{locs}{locData from a GenoSet object}
}
\value{
  DataFrame of Rle objects with nrows same as locs and one
  column for each sample
}
\description{
  Given segments, make a DataFrame of Rle objects for each
  sample
}
\details{
  Take table of segments from CBS, convert DataTable of Rle
  objects for each sample.
}
\examples{
data(genoset)
  seg.list = runCBS( lrr(baf.ds), locData(baf.ds), return.segs=TRUE )
  segs2RleDataFrame( seg.list, locData(baf.ds) )  # Loop segs2Rle on list of data.frames in seg.list
}
\author{
  Peter Haverty
}
\seealso{
  Other "segmented data": \code{\link{bounds2Rle}},
  \code{\link{rangeSegMeanLength}},
  \code{\link{rangeSegMeanLength}},
  \code{\link{rangeSegMeanLength}}, \code{\link{runCBS}},
  \code{\link{segPairTable}}, \code{\link{segPairTable}},
  \code{\link{segPairTable}}, \code{\link{segs2Granges}},
  \code{\link{segs2RangedData}}, \code{\link{segs2Rle}},
  \code{\link{segTable}}, \code{\link{segTable}},
  \code{\link{segTable}}

  Other segments: \code{\link{segs2RangedData}}
}

