\docType{methods}
\name{genoPlot}
\alias{genoPlot}
\alias{genoPlot-methods}
\alias{genoPlot,numeric,numeric-method}
\alias{genoPlot,numeric,Rle-method}
\alias{genoPlot,RangedDataOrGenoSetOrGenomicRanges,ANY-method}
\title{Plot data along the genome}
\arguments{
  \item{x}{GenoSet (or descendant), RangedData, or GRanges}

  \item{y}{numeric or Rle}

  \item{locs}{RangedData, like locData slot of GenoSet}

  \item{chr}{Chromosome to plot, NULL by default for full
  genome}

  \item{add}{Add plot to existing plot}

  \item{xlab}{character, label for x-axis of plot}

  \item{ylab}{character, label for y-axis of plot}

  \item{col}{character, color to plot lines or points}

  \item{lwd}{numeric, line width for segment plots from an
  Rle}

  \item{pch}{character or numeric, printing character, see
  points}

  \item{xlim}{integer, length two, bounds for genome
  positions. Used in conjunction with "chr" to subset data
  for plotting.}

  \item{...}{Additional plotting args}
}
\value{
  nothing
}
\description{
  Plot location data and chromosome boundaries from a
  GenoSet, RangedData, or GRanges object against data from
  a numeric or Rle. Specifying a chromosome name and
  optionally a 'xlim' will zoom into one chromosome region.
  If more than one chromosome is present, the chromosome
  boundaries will be marked. Alternatively, for a numeric x
  and a numeric or Rle y, data in y can be plotted at
  genome positions x. In this case, chromosome boundaries
  can be taken from the argument locs. If data for y-axis
  comes from a Rle lines are plotted representing segments.
  X-axis tickmarks will be labeled with genome positions in
  the most appropriate units.
}
\section{Methods}{
  \describe{

  \item{\code{signature(x =
  "RangedDataOrGenoSetOrGenomicRanges", y = "ANY")}}{ Plot
  feature locations and data from one sample. }

  \item{\code{signature(x = "numeric", y = "numeric")}}{
  Plot numeric location and a vector of numeric data. }

  \item{\code{signature(x = "numeric", y = "Rle")}}{ Plot
  numeric location and a vector of Rle data. Uses lines for
  Rle runs. } }
}
\examples{
data(genoset)
genoPlot( x=genoset.ds,y=genoset.ds[,1,"lrr"] )
genoPlot( genoPos(genoset.ds), genoset.ds[,1,"lrr"], locs=locData(genoset.ds) ) # The same
genoPlot( 1:10, Rle(c(rep(0,5),rep(3,4),rep(1,1))) )
}
\author{
  Peter M. Haverty
}
\seealso{
  Other "genome plots": \code{\link{genomeAxis}}
}

