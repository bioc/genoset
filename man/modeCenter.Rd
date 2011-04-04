\name{modeCenter}
\alias{modeCenter}
\title{Center continuous data on mode...}
\usage{modeCenter(ds)
}
\description{Center continuous data on mode}
\details{Copynumber data distributions are generally multi-modal. It is often assumed that
the tallest peak represents "normal" and should therefore be centered on a
log2ratio of zero. This function uses the density function to find the mode of
the dominant peak and subtracts that value from the input data.}
\value{numeric matrix}
\author{Peter M. Haverty}
\arguments{\item{ds}{numeric matrix}
}
\examples{modeCenter( matrix( rnorm(150, mean=0), ncol=3 ))}

