\name{segTable}
\alias{segTable}
\title{Take a DataFrame of Rle vectors and make a list of data...}
\usage{segTable(df, locs)
}
\description{Take a DataFrame of Rle vectors and make a list of data.frames}
\details{Like the inverse of segs2RleDataFrame. Take a DataFrame with Rle
columns and the locData RangedData both from a GenoSet object
and make a list of data.frames each like the result of CBS's
segment.  Note the loc.start and loc.stop will correspond
exactly to probe locations in locData and the input to
segs2RleDataFrame are not necessarily so.}
\value{list of data.frames with columns ID, chrom, loc.start, loc.end, num.mark, seg.mean}
\author{Peter M. Haverty}
\arguments{\item{df}{list or DataFrame of Rle vectors}
\item{locs}{RangedData with rows corresponding to rows of df}
}
\examples{data(genoset)
seg.list = runCBS( lrr(baf.ds), locData(baf.ds), return.segs=TRUE )
df = segs2RleDataFrame( seg.list, locData(baf.ds) )  # Loop segs2Rle on list of data.frames in seg.list
assayDataElement( baf.ds, "lrr.segs" ) = df
segTable( df, locData(baf.ds) )
segTable( assayDataElement(baf.ds,"lrr.segs"), locData(baf.ds) )}

