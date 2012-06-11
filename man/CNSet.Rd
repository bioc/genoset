\name{CNSet}
\alias{CNSet}
\title{Create a CNSet object}
\usage{
  CNSet(locData, cn = NULL, pData = NULL, annotation = "",
    assayData = NULL, ...)
}
\arguments{
  \item{locData}{A GRanges or RangedData object specifying
  feature chromosome locations. featureNames (names or
  rownames) are required to match featureNames of
  matrices.}

  \item{cn}{numeric matrix of copy number data with
  rownames matching featureNames and colnames matching
  sampleNames}

  \item{pData}{A data frame with rownames matching all data
  matrices}

  \item{annotation}{character, string to specify
  chip/platform type}

  \item{assayData}{assayData, usually an environment}

  \item{...}{More matrix or DataFrame objects to include in
  assayData}
}
\value{
  A CNSet object
}
\description{
  This function is the preferred method for creating a new
  CNSet object. Users are generally discouraged from
  calling "new" directly. This CNSet function enforces the
  requirement for a "cn" matrix.  This and any other "..."
  arguments will become part of the assayData slot of the
  resulting object. "..." can be matrices or DataFrame
  objects (from the IRanges package). This function passes
  control to the "initGenoSet" method which performs
  argument checking including dimname matching among
  relevant slots and sets everything to genome order.
  Genome order can be disrupted by "[" or "[[" calls and
  will be checked by methods that require it.
}
\examples{
test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]
joe = CNSet(
   locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
   cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
   pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
   annotation="SNP6"
   )
}
\author{
  Peter M. Haverty
}

