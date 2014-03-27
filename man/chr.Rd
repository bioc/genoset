\name{chr}
\alias{chr}
\title{Look up chromosome for each feature}
\arguments{
  \item{object}{GRanges, RangedData or GenoSet}
}
\value{
  character vector of chromosome positions for each feature
}
\description{
  Chromosome name for each feature
}
\details{
  Get chromosome name for each feature.  Returns character,
  not the factor 'space'.
}
\examples{
data(genoset)
  chr(genoset.ds)  # c("chr1","chr1","chr1","chr1","chr3","chr3","chrX","chrX","chrX","chrX")
  chr(locData(genoset.ds))  # The same
}

