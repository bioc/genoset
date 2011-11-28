\name{segs2Rle}
\alias{segs2Rle}
\title{Make Rle from segments for one sample...}
\usage{segs2Rle(segs, locs)}
\description{Make Rle from segments for one sample}
\details{Take output of CBS, make Rle representing all features in 'locs' ranges. CBS output contains
run length and run values for genomic segmetns, which could very directly be converted into a Rle.
However, as NA values are often removed, especially for mBAF data, these run lengths do not
necessarily cover all features in every sample. Using the start and top positions of each segment 
and the location of each feature, we can make a Rle that represents all features.}
\value{Rle with run lengths and run values covering all features in the data set.}
\author{Peter M. Haverty \email{phaverty@gene.com}}
\arguments{\item{segs}{data.frame of segments, formatted as output of segment function from DNAcopy package}
\item{locs}{RangedData, like locData slot of a GenoSet}}
\examples{data(genoset)
segs = runCBS( lrr(baf.ds), locData(baf.ds), return.segs=TRUE )
segs2Rle( segs[[1]], locData(baf.ds) )  # Take a data.frame of segments, say from DNAcopy's segment function, and make Rle's using probe locations in the RangedData locs}
