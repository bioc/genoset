\name{chrNames}
\alias{chrNames}
\alias{chrNames,GenoSet-method}
\alias{chrNames,GRanges-method}
\alias{chrNames,RangedData-method}
\alias{chrNames<-}
\alias{chrNames<-,GenoSet-method}
\alias{chrNames<-,GRanges-method}
\alias{chrNames<-,RangedData-method}
\title{Get list of unique chromosome names}
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
\author{
  Peter M. Haverty
}

