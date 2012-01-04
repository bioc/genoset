\docType{package}
\name{genoset-package}
\alias{genoset}
\alias{genoset-package}
\title{GenoSet: An eSet for data with genome locations}
\description{
  Load, manipulate, and plot copynumber and BAF data.
  GenoSet class extends eSet by adding a "locData" slot for
  a RangedData object from the IRanges package. This object
  contains feature genome location data and provides for
  simple subsetting on genome location. CNSet and BAFSet
  extend GenoSet and require assayData matrices for Copy
  Number (cn) or Log-R Ratio (lrr) and B-Allele Frequency
  (baf) data. Implements and provides convenience functions
  for processing of copy number and B-Allele Frequency
  data.
}

