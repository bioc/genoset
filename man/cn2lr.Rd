\name{cn2lr}
\alias{cn2lr}
\alias{cn2lr-methods}
\title{Take vector or matrix of copynumber values, convert to log2ratios
Utility function for converting copynumber units (2 is normal) to log2ratio units (two is normal). If ploidy
is provided lr is log2(cn/ploidy), otherwise log2(cn/2).}
\usage{
cn2lr(x, ploidy)
}
\arguments{
  \item{x}{numeric vector or matrix, or DataFrame with
  numeric-like columns (Rle typicaly). Assumed to be in
  copynumber units.}

  \item{ploidy}{numeric, of length ncol(x). Ploidy of each
  sample.}
}
\value{
data of same type as "x" transformed into log2ratio units
}
\description{
Take vector or matrix of copynumber values, convert to
log2ratios Utility function for converting copynumber units
(2 is normal) to log2ratio units (two is normal). If ploidy
is provided lr is log2(cn/ploidy), otherwise log2(cn/2).
}
\author{
Peter M. Haverty \email{phaverty@gene.com}
}
\seealso{
lr2cn
}

