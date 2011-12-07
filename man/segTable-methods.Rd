\docType{methods}
\name{segTable}
\alias{segTable}
\alias{segTable,DataFrame-method}
\alias{segTable,Rle-method}
\title{Convert Rle objects to tables of segments}
\arguments{
  \item{object}{Rle or list/DataFrame of Rle vectors}

  \item{locs}{RangedData with rows corresponding to rows of
  df}

  \item{sample.name}{character for single Rle optionally
  include "ID" column with this sample name}
}
\value{
  one or a list of data.frames with columns ID, chrom,
  loc.start, loc.end, num.mark, seg.mean
}
\description{
  Like the inverse of segs2Rle and segs2RleDataFrame. Takes
  a Rle or a DataFrame with Rle columns and the locData
  RangedData both from a GenoSet object and make a list of
  data.frames each like the result of CBS's segment.  Note
  the loc.start and loc.stop will correspond exactly to
  probe locations in locData and the input to
  segs2RleDataFrame are not necessarily so.
}
\examples{
data(genoset)
  seg.list = runCBS( lrr(baf.ds), locData(baf.ds), return.segs=TRUE )
  df = segs2RleDataFrame( seg.list, locData(baf.ds) )  # Loop segs2Rle on list of data.frames in seg.list
  assayDataElement( baf.ds, "lrr.segs" ) = df
  segTable( df, locData(baf.ds) )
  segTable( assayDataElement(baf.ds,"lrr.segs"), locData(baf.ds) )
  segTable( assayDataElement(baf.ds,"lrr.segs")[,1], locData(baf.ds), sampleNames(baf.ds)[1] )
}
\author{
  Peter M. Haverty
}
\seealso{
  Other "segmented data": \code{\link{runCBS}},
  \code{\link{segs2RangedData}}, \code{\link{segs2Rle}},
  \code{\link{segs2RleDataFrame}}
}

