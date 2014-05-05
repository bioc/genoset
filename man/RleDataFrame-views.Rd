\name{.do_rledf_views}
\alias{.do_rledf_views}
\title{Calculate min/max/sum/mean/whichmin/whichmax over each view on each column of an RleDataFrame.}
\usage{
.do_rledf_views(x, ir, na.rm = FALSE, simplify = TRUE, FUN,
  FUN.TYPE = c("numeric", "double", "integer", "logical"))
}
\arguments{
  \item{x}{RleDataFrame}

  \item{ir}{IRanges or matrix, views on every Rle. If ir is
  a matrix, it is converted to an IRanges using the first
  two columns as the starts and stops. Names for the
  IRanges are taken from the rownames of the matrix.}

  \item{na.rm}{scalar logical, ignore NAs in calculations?}

  \item{simplify}{scalar logical, simplify result? For
  single view, vector, otherwise matrix with one row per
  view.}

  \item{FUN}{S4Generic or scalar character with the name of
  an S4 Generic}

  \item{FUN.TYPE}{scalar character, the storage mode for
  the returned vector or matrix (when simplify==TRUE).}
}
\value{
With simplify == TRUE, vector for single view or single
column x with simplify == TRUE, matix otherwise. When
simplify == FALSE, a list of length ncol(x).
}
\description{
Loop over Rles in RleDataFrame, calculate the appropriate
statistic for each view. If simplify == FALSE, returns a
vector for each Rle. If simplify == TRUE, returns a vector
for the case of a single view, otherwise, a matrix.
Rownames for the matrix are taken from the names of the
argument \code{ir}.
}
\keyword{internal}

