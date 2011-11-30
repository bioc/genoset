\name{BAFSet}
\alias{BAFSet}
\title{Create a BAFSet object}
\usage{
  BAFSet(locData, lrr, baf, pData = NULL, annotation = "",
    universe = NULL, ...)
}
\arguments{
  \item{locData}{A RangedData object specifying feature
  chromosome locations. Rownames are required to match
  featureNames.}

  \item{lrr}{numeric matrix of copy number data with
  rownames matching sampleNames and colnames matching
  sampleNames}

  \item{baf}{numeric matrix of B-Allele Frequency data with
  rownames matching sampleNames and colnames matching
  sampleNames}

  \item{pData}{A data frame with rownames matching all data
  matrices}

  \item{annotation}{character, string to specify
  chip/platform type}

  \item{universe}{character, a string to specify the genome
  universe for locData}

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
\examples{
test.sample.names = LETTERS[11:13]
  probe.names = letters[1:10]
  locData.rd = RangedData(ranges=IRanges(start=c(1,4,3,2,5:10),width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18")
  bs = BAFSet(
    locData=locData.rd,
    lrr=matrix(1:30,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    baf=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
    pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
    annotation="SNP6"
)
}
\author{
  Peter M. Haverty
}
\seealso{
  bafset-class, genoset-class
}

