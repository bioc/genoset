\name{chrNames}
\alias{chrNames}
\title{Get list of unique chromosome names}
\usage{
chrNames(object)
}
\arguments{
  \item{object}{RangedData or GenoSet}
}
\value{
character vector with names of chromosomes
}
\description{
Get list of unique chromosome names
}
\examples{
data(genoset)
  chrNames(genoset.ds) # c("chr1","chr3","chrX")
  chrNames(locData(genoset.ds))  # The same
  chrNames(genoset.ds) = sub("^chr","",chrNames(genoset.ds))
}

