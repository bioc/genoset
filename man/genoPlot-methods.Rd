\docType{methods}
\name{genoPlot}
\alias{genoPlot}
\alias{genoPlot-methods}
\alias{genoPlot,BAFSet,ANY-method}
\alias{genoPlot,CNSet,ANY-method}
\alias{genoPlot,GenoSet,ANY-method}
\alias{genoPlot,numeric,numeric-method}
\alias{genoPlot,numeric,Rle-method}
\title{Plot data along the genome}
\arguments{
  \item{x}{GenoSet (or descendant) or numeric with
  chromosome or genome positions. See Details.}

  \item{y}{numeric or Rle, or single integer. See Details.}

  \item{element}{character, when x is a GenoSet, the y-th
  column of this assayDataElement is used for the y-axis
  data. The "cn" element is the default for a CNSet.  For a
  BAFSet, it is "lrr".}

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
\details{
  genoPlot has two modes.  When x and y are simple data,
  data in y is plotted against locations in x. If location
  information is provided in the "locs" argument,
  chromosome bounds and genome scale x-axis tick labels
  will be added.  When x is a Genoset-derived object,
  genome location info is extracted for the x-values and
  the y values come from the y-th column in the
  assayDataElement listed in the "element" argument (i.e.
  x[ , y, element]).
}
\section{Methods}{
  \describe{

  \item{\code{signature(x = "GenoSet", y = "ANY")}}{ Plot
  feature locations and data from one sample. }

  \item{\code{signature(x = "BAFSet", y = "ANY")}}{ Plot
  feature locations and data from one sample. Defaults to
  data from 'lrr' assayDataElement. }

  \item{\code{signature(x = "CNSet", y = "ANY")}}{ Plot
  feature locations and data from one sample. Defaults to
  data from 'cn' assayDataElement. }

  \item{\code{signature(x = "numeric", y = "numeric")}}{
  Plot numeric location and a vector of numeric data. }

  \item{\code{signature(x = "numeric", y = "Rle")}}{ Plot
  numeric location and a vector of Rle data. Uses lines for
  Rle runs. } }
}
\examples{
data(genoset)
genoPlot( baf.ds,1,element="lrr")
genoPlot( genoPos(baf.ds), assayDataElement(baf.ds,"lrr")[,1], locs=locData(baf.ds) ) # The same
genoPlot( 1:10, Rle(c(rep(0,5),rep(3,4),rep(1,1))) )
}
\author{
  Peter M. Haverty
}
\seealso{
  Other "genome plots": \code{\link{genomeAxis}}
}

