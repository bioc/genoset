\name{gcCorrect}
\alias{gcCorrect}
\title{cgCorrect}
\usage{gcCorrect(ds, gc, retain.mean=TRUE)
}
\description{Correct copy number for GC content}
\details{Copy number estimates from various platforms show "Genomic Waves" (Diskin et al.,
Nucleic Acids Research, 2008) where copy number trends with local GC content.
This function regresses copy number on GC percentage and removes the effect
(returns residuals). GC content should be smoothed along the genome in wide
windows >= 100kb.}
\value{numeric matrix, residuals of ds regressed on gc}
\author{Peter M. Haverty}
\arguments{\item{ds}{numeric matrix of copynumber or log2ratio values, samples in columns}
\item{gc}{numeric vector, GC percentage for each row of ds, must not have NAs}
\item{retain.mean}{logical, center on zero or keep same mean?}
}
\examples{gc = runif(n=100, min=1, max=100)
ds = rnorm(100) + (0.1 * gc)
gcCorrect(ds, gc)}

