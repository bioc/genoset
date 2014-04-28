\docType{methods}
\name{colMeans}
\alias{colMeans}
\alias{colMeans,DataFrame-method}
\alias{colMeans,RleDataFrame-method}
\alias{colSums}
\alias{colSums,RleDataFrame-method}
\alias{rowMeans}
\alias{rowMeans,ANY-method}
\alias{rowMeans,RleDataFrame-method}
\title{Means of columns}
\usage{
colMeans(x, na.rm = FALSE, dims = 1, ...)

\S4method{colMeans}{RleDataFrame}(x, na.rm = TRUE)

\S4method{colMeans}{DataFrame}(x, na.rm = TRUE, dims = 1L)

colSums(x, na.rm = FALSE, dims = 1, ...)

\S4method{colSums}{RleDataFrame}(x, na.rm = TRUE)

rowMeans(x, na.rm = FALSE, dims = 1, ...)

\S4method{rowMeans}{ANY}(x, na.rm = FALSE, dims = 1L)

\S4method{rowMeans}{RleDataFrame}(x, na.rm = FALSE, dims = 1L)
}
\arguments{
  \item{x}{DataFrame}

  \item{na.rm}{logical}

  \item{dims}{integer}

  \item{...}{in generic, for extra args in methods.}
}
\description{
Calculate means of columns of a DataFrame as if it were a
matrix. Allows colmeans in rangeSampleMeans for DataTable
just like a real matrix. I'm sure there is much more clever
way to do this using aggregate.
}
\examples{
df.ds = DataFrame( a = Rle(c(5,4,3),c(2,2,2)), b = Rle(c(3,6,9),c(1,1,4)) )
 mat.ds = matrix( c(5,5,4,4,3,3,3,6,9,9,9,9), ncol=2, dimnames=list(NULL,c("a","b")))
 \dontrun{ identical( colMeans(df.ds), colMeans(mat.ds) ) }
}

