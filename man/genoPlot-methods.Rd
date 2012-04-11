\docType{methods}
\name{genoPlot}
\alias{genoPlot}
\alias{genoPlot,BAFSet,ANY-method}
\alias{genoPlot,CNSet,ANY-method}
\alias{genoPlot,GenoSet,ANY-method}
\alias{genoPlot,numeric,numeric-method}
\alias{genoPlot,numeric,Rle-method}
\title{Plot data along the genome}
\arguments{
  \item{sample}{A index or sampleName to plot}

  \item{element}{character, name of element in assayData to
  plot}

  \item{x}{GenoSet (or descendant) or numeric with
  chromosome or genome positions}

  \item{y}{numeric or Rle, values to be used for
  y-dimension, run start and stop indices or numeric with
  all values mapped to values in x for x-dimension or index
  of sample to be plotted if x is a GenoSet.}

  \item{element}{character, when x is a GenoSet, the name
  of the assayDataElement to plot from.}

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

  \item{...}{Additional plotting args}
}
\value{
  nothing
}
\description{
  For a GenoSet object, data for a specified sample in a
  specified assayDataElement can be plotted along the
  genome.  One chromosome can be specified if desired. If
  more than one chromosome is present, the chromosome
  boundaries will be marked. Alternatively, for a numeric x
  and a numeric or Rle y, data in y can be plotted at
  genome positions y. In this case, chromosome boundaries
  can be taken from the argument locs. If data for y-axis
  comes from a Rle, either specified directly or coming
  from the specified assayData element and sample, lines
  are plotted representing segments.
}
\examples{
data(genoset)
  genoPlot( baf.ds,1,element="lrr")
  genoPlot( genoPos(baf.ds), assayDataElement(baf.ds,"lrr")[,1], locs=locData(baf.ds) ) # The same
  genoPlot( 1:10, Rle(c(rep(0,5),rep(3,4),rep(1,1))) )
}
\author{
  Peter M. Haverty

  Peter M. Haverty
}
\seealso{
  Other "genome plots": \code{\link{genomeAxis}}
}

