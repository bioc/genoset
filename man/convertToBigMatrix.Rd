\name{convertToBigMatrix}
\alias{convertToBigMatrix}
\title{Make standard matrices in a GenoSet filebacked bigmatrix objects}
\usage{
  convertToBigMatrix(object, prefix = "bigmat",
    path = "bigmat")
}
\arguments{
  \item{object}{GenoSet}

  \item{prefix}{character, prefix for all bigmatrix related
  files}

  \item{path}{character, directory to be created for all
  bigmatrix files, can be pre-existing.}
}
\value{
  GenoSet or related, updated copy of "object"
}
\description{
  Make standard matrices in a GenoSet filebacked bigmatrix
  objects. Something like a factor can be obtained using
  integer assayDataElements with a "levels" attribute. The
  levels attribute will be maintained. Such objects will be
  stored as char on disk if there are < 128 levels, and
  integer otherwise. "nlevels" and "levels" will work on
  these objects as they only require the levels attribute.
  The "as.character" functionality of a factor can be
  obtained like this: levels(assayDataElement(ds,"geno"))[
  ds[1:5,1:5,"geno"] ] for a GenoSet called "ds" with a
  factor-like element called "geno".
}
\examples{
\dontrun{ ds = convertToBigMatrix(ds) }
}
\author{
  Peter M. Haverty \email{phaverty@gene.com}
}

