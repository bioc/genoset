\name{colMeans}
\alias{colMeans-methods}
\alias{colMeans}
\title{Means of columns...}
\description{Means of columns}
\details{\code{colMeans-methods}: Get means of columns of a DataFrame as if it were a matrix

}
\author{Peter M. Haverty}
\arguments{\item{x}{DataFrame}
\item{na.rm}{logical}
\item{dims}{integer}
}
\examples{df.ds = DataFrame( a = Rle(c(5,4,3),c(2,2,2)), b = Rle(c(3,6,9),c(1,1,4)) )
mat.ds = matrix( c(5,5,4,4,3,3,3,6,9,9,9,9), ncol=2, dimnames=list(NULL,c("a","b")))
identical( colMeans(df.ds), colMeans(mat.ds) )}
\alias{colMeans,-method}
\alias{colMeans}

