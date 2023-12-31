% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segments.R
\docType{methods}
\name{segTable}
\alias{segTable}
\alias{segTable,Rle-method}
\alias{segTable,DataFrame-method}
\title{Convert Rle objects to tables of segments}
\usage{
segTable(object, ...)

\S4method{segTable}{Rle}(object, locs = NULL, chr.ind = NULL,
  start = NULL, end = NULL, factor.chr = TRUE)

\S4method{segTable}{DataFrame}(object, locs, factor.chr = TRUE,
  stack = FALSE)
}
\arguments{
\item{object}{Rle or RleDataFrame}

\item{...}{in generic, for extra args in methods}

\item{locs}{GenomicRanges with rows corresponding to rows of df}

\item{chr.ind}{matrix, like from chrIndices method}

\item{start}{integer, vector of feature start positions}

\item{end}{integer, vector of feature end positions}

\item{factor.chr}{scalar logical, make 'chrom' column a factor?}

\item{stack}{logical, rbind list of segment tables for each sample and add 'Sample' column?}
}
\value{
one or a list of data.frames with columns chrom, loc.start, loc.end, num.mark, seg.mean
}
\description{
Like the inverse of segs2Rle and segs2RleDataFrame. Takes a
Rle or a RleDataFrame and the rowRanges both from a GenoSet object
and makes a list of data.frames each like the result of CBS's
segment.  Note the loc.start and loc.stop will correspond
exactly to probe locations in rowRanges and the input to
segs2RleDataFrame are not necessarily so. For a DataFrame, the
argument \code{stack} combines all of the individual data.frames
into one large data.frame and adds a 'Sample' column of sample ids.
}
\details{
For a Rle, the user can provide \code{locs} or \code{chr.ind},
\code{start} and \code{stop}.  The latter is surprisingly much faster
and this is used in the DataFrame version.
}
\examples{
  data(genoset,package='genoset')
  seg.list = runCBS( genoset.ds[, , 'lrr'], rowRanges(genoset.ds), return.segs=TRUE )
  df = segs2RleDataFrame( seg.list, rowRanges(genoset.ds) )  # Loop segs2Rle on list of data.frames in seg.list
  genoset.ds[ , , 'lrr.segs'] = df
  segTable( df, rowRanges(genoset.ds) )
  segTable( genoset.ds[ , , 'lrr.segs'], rowRanges(genoset.ds) )
  segTable( genoset.ds[ , 1, 'lrr.segs'], rowRanges(genoset.ds), colnames(genoset.ds)[1] )
}
\seealso{
Other 'segmented data': \code{\link{bounds2Rle}},
  \code{\link{rangeSegMeanLength}}, \code{\link{runCBS}},
  \code{\link{segPairTable}}, \code{\link{segs2Granges}},
  \code{\link{segs2RleDataFrame}}, \code{\link{segs2Rle}}
}
\concept{'segmented data'}
