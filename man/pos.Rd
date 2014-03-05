\name{pos}
\alias{pos}
\alias{pos,RangedDataOrGenoSetOrGenomicRanges-method}
\title{Positions for features}
\usage{
pos(object)
}
\arguments{
  \item{object}{GRanges, RangedData or GenoSet}
}
\value{
numeric vector of feature positions within a chromosome
}
\description{
Chromosome position of features
}
\details{
Get chromosome position of features/ranges. Defined as
floor of mean of start and end.
}
\examples{
data(genoset)
  pos(genoset.ds)  # 1:10
  pos(locData(genoset.ds))  # The same
}
\author{
Peter Haverty
}

