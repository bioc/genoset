\name{RleDataFrame-views}
\alias{RleDataFrame-views}

\alias{rangeSums,RleDataFrame-method}
\alias{rangeMeans,RleDataFrame-method}
\alias{rangeMins,RleDataFrame-method}
\alias{rangeMaxs,RleDataFrame-method}
\alias{rangeWhichMins,RleDataFrame-method}
\alias{rangeWhichMaxs,RleDataFrame-method}

\title{Calculate summary statistics on views of an RleDataFrame}

\description{
  
}

\usage{
rangeSums(x, ir, na.rm=TRUE, simplify=TRUE)
rangeMeans(x, ir, na.rm=TRUE, simplify=TRUE)
rangeMins(x, ir, na.rm=TRUE, simplify=TRUE)
rangeMaxs(x, ir, na.rm=TRUE, simplify=TRUE)
rangeWhichMins(x, ir, na.rm=TRUE, simplify=TRUE)
rangeWhichMaxs(x, ir, na.rm=TRUE, simplify=TRUE)
}

\arguments{
  \item{x}{
    RleDataFrame
  }
  \item{ir}{
    Matrix with two columns or IRanges representing ranges of rows of
  \code{x} to process. if \code{ir} is a matrix, an IRanges is
  constructed assuming the first two columns represent the start and end
  of the ranges. The names for the IRanges is taken from the rownames of
  the matrix. Such a matrix can constructed with \code{boundingIndicesByChr}.
  }
  \item{na.rm}{
    Scalar logical. Ignore NAs in calculations?
  }
  \item{simplify}{
    Scalar logical. Simplify result? If TRUE, the return value will be a
  vector or matrix. For a single view, a vector will be
  returned. Otherwise a matrix with one row per view and one column per
  column of x will be returned. If FALSE, the return value will be a
  list of length \code{ncol(x)} of vectors of length \code{nrow(ir)}.
  }
}

\value{
   With \code{simplify == TRUE}, a vector for single view or a matrix
otherwise. When \code{simplify == FALSE}, a list of vectors length ncol(x) where each element is of length \code{nrows(ir)}.
}

\examples{

df = RleDataFrame(list(a=Rle(1:5, rep(2, 5))), b=Rle(1:5, rep(2, 5)),
row.names=LETTERS[1:10])
mat = matrix(c(1,4,3,5),ncol=2,dimnames=list(c("Gene1","Gene2"),c("start","end")))
ir = IRanges(start=c(1, 4), end=c(3, 5), names=c("Gene1","Gene2"))

rangeMeans(df,ir,simplify=FALSE)
rangeMeans(df,ir,simplify=TRUE)
rangeMeans(df,mat,simplify=TRUE)

rangeSums(df,ir,simplify=TRUE)
rangeMins(df,ir,simplify=TRUE)
rangeMaxs(df,ir,simplify=TRUE)
rangeWhichMins(df,ir,simplify=TRUE)
rangeWhichMaxs(df,ir,simplify=TRUE)
  
}

\seealso{
  \linkS4class{RleDataFrame} \link{boundingIndicesByChr}
}

\keyword{methods}
