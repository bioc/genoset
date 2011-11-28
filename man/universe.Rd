\name{universe}
\alias{universe}
\title{Get and set the genome universe annotation.}
\description{Genome universe for locData}
\details{The genome positions of the features in locData. The UCSC notation (e.g. hg18, hg19, etc.) should be used.}
\value{character, e.g. hg19}
\author{Peter M. Haverty}
\arguments{\item{x}{GenoSet}}
\examples{data(genoset)
universe(genoset.ds)
universe(genoset.ds) = "hg19"}
