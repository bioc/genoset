\name{chrOrder}
\alias{chrOrder}
\title{Order chromosome names in proper genome order...}
\usage{chrOrder(chr.names)
}
\description{Order chromosome names in proper genome order}
\details{Chromosomes make the most sense orded by number, then by letter.}
\value{character vector of chromosome names in proper order}
\author{Peter M. Haverty}
\arguments{\item{chr.names}{character, vector of unique chromosome names}
}
\examples{chrOrder(c("chr5","chrX","chr3","chr7","chrY"))  #  c("chr3","chr5","chr7","chrX","chrY")}

