\name{.do_rledf_views}
\alias{.do_rledf_views}
\title{Calculate min/max/sum/mean/whichmin/whichmax over each view on each column of an RleDataFrame.}
\usage{
.do_rledf_views(x, ir, na.rm = FALSE, simplify = TRUE, FUN,
  FUN.TYPE = c("numeric", "double", "integer", "logical"))
}
\arguments{
  \item{x}{RleDataFrame}

  \item{ir}{IRanges or matrix, views on every Rle. If
  \code{ir} is a matrix, it is converted to an IRanges
  using the first two columns as the starts and stops.
  Names for the IRanges are taken from the rownames of the
  matrix. Such a matrix can be constructed with
  \code{boundingIndicesByChr}.}

  \item{na.rm}{scalar logical, ignore NAs in calculations?}

  \item{simplify}{scalar logical, simplify result? For
  single view or Rle, a vector, otherwise a matrix with one
  row per view.}

  \item{FUN}{S4 Generic function or scalar character with
  the name of an S4 Generic function}

  \item{FUN.TYPE}{scalar character, the storage mode for
  the returned vector or matrix (when simplify==TRUE).}
}
\value{
With \code{simplify == TRUE}, a vector for single view or
single Rle or a matrix otherwise. When \code{simplify ==
FALSE}, a list of vectors length ncol(x) where each element
is of length \code{nrows(ir)}.
}
\description{
Loop over the Rle objects in an RleDataFrame, calculate the
appropriate statistic for each view. If simplify == FALSE,
this function returns a vector for each Rle. If simplify ==
TRUE, it returns a vector for the case of a single view or
Rle, otherwise, a matrix. Rownames for the matrix are taken
from the names of the argument \code{ir}.
}
\seealso{
RleDataFrame boundingIndicesByChr
}
\keyword{internal}

