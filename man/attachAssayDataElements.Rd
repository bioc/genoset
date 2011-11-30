\name{attachAssayDataElements}
\alias{attachAssayDataElements}
\title{Attach on-disk matrices into assayData}
\usage{
  attachAssayDataElements(object)
}
\arguments{
  \item{object}{eSet}
}
\value{
  assayData in storage mode of input assayData, invisibly.
  Re-assignment back original eSet only necessary if using
  a list type assayData.
}
\description{
  GenoSet objects can hold big.matrix objects in their
  assayData slot environment. After re-loading the GenoSet
  from disk, these objects will each need to be re-attached
  to their on-disk component using their resource locators
  stored in their "desc" attributes. This function checks
  each assayDataElement to see if it is an un-attached
  big.matrix object, re-attaching if necessary. All other
  assayDataElements are left untouched. In later releases
  this function will also handle other on-disk types, like
  HDF5-based matrices.
}
\details{
  *** Intentional side-effects *** Environment type
  assayData objects, even "lockedEnvironment" objects, will
  be updated in place (same pointer). This allows for
  functions trying to access assayDataElements to attach
  before access, rather than crashing R.
}
\author{
  Peter M. Haverty \email{phaverty@gene.com}
}

