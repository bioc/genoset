\name{genoset-deprecated}
\Rdversion{1.1}
\alias{genoset-deprecated}

\title{Functions deprecated in genoset}
\description{
  Several functions have been deprecated in an effort to allow either a
  GRanges or a RangedData in the locData slot. RangedData-specific
  functions have been replaced with more generically named functions.
}
\section{Deprecated Functions}{
  \describe{
    \item{space}{space(locData(object)) or seqnames(locData(object)) will get you what you need for a RangedData or GRanges, respectively.}
    \item{ranges}{Work on locData(object) directly with ranges(locData(object)).}
    \item{names}{Please try chrNames.}
  }
}
\author{
  Peter M. Haverty \email{phaverty@gene.com}
}
\keyword{deprecated}
