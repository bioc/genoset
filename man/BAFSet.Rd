\name{BAFSet}
\alias{BAFSet}
\alias{BAFSet-defunct}
\title{Create a BAFSet object}
\usage{
  BAFSet(locData, lrr = NULL, baf = NULL, pData = NULL,
    annotation = "", universe, assayData = NULL, ...)
}
\arguments{
  \item{locData}{A GRanges or RangedData object specifying
  feature chromosome locations. featureNames (names or
  rownames) are required to match featureNames of
  assayData.}

  \item{lrr}{numeric matrix of copy number data with
  rownames matching featureNames and colnames matching
  sampleNames}

  \item{baf}{numeric matrix of B-Allele Frequency data with
  rownames matching featureNames and colnames matching
  sampleNames}

  \item{pData}{A data frame with rownames matching all data
  matrices}

  \item{annotation}{character, string to specify
  chip/platform type}

  \item{universe}{character, a string to specify the genome
  universe for locData. Overrides any universe/genome data
  in locData.}

  \item{assayData}{assayData, usually an environment}

  \item{...}{More matrix or DataFrame objects to include in
  assayData slot}
}
\value{
  A BAFSet object
}
\description{
  This function is the preferred method for creating a new
  BAFSet object. Users are generally discouraged from
  calling "new" directly. This BAFSet function enforces the
  requirement for "lrr" and "baf" matrices.  These and any
  other "..." arguments will become part of the assayData
  slot of the resulting object. "..." can be matrices or
  DataFrame objects (from the IRanges package). This
  function passes control to the "initGenoSet" method which
  performs argument checking including dimname matching
  among relevant slots and sets everything to genome order.
  Genome order can be disrupted by "[" or "[[" calls and
  will be checked by methods that require it.
}
\details{
  The BAFSet class is defunct. Please use GenoSet. BAFSet
  only added the baf/lrr getter/setter functions, which are
  redundant with x[, , 'baf'] and x[, , 'lrr'] now.
}
\author{
  Peter M. Haverty
}
\seealso{
  bafset-class, genoset-class
}

