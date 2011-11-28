\name{space}
\alias{space}
\title{Get space factor for GenoSet...}
\description{Get space factor for GenoSet}
\details{locData slot holds a RangedData, which keeps the chromsome of each
feature in a factor names 'space'.}
\value{factor}
\author{Peter M. Haverty}
\arguments{\item{x}{GenoSet}}
\examples{data(genoset)
space(genoset.ds)
start(genoset.ds)
end(genoset.ds)
names(genoset.ds)
ranges(genoset.ds) # Returns a RangesList
elementLengths(genoset.ds) # Returns the number of probes per chromosome}
