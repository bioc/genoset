\name{rangeSampleMeans}
\alias{rangeSampleMeans}
\title{Average features in ranges per sample}
\usage{
  rangeSampleMeans(query.rd, subject, assay.element)
}
\arguments{
  \item{query.rd}{RangedData object representing genomic
  regions (genes) to be averaged.}

  \item{subject}{A GenoSet object or derivative}

  \item{assay.element}{character, name of element in
  assayData to use to extract data}
}
\value{
  numeric matrix of features in each range averaged by
  sample
}
\description{
  This function takes per-feature genomic data and returns
  averages for each of a set of genomic ranges. The most
  obvious application is determining the copy number of a
  set of genes. The features corresponding to each gene are
  determined with boundingIndices such that all features
  with the bounds of a gene (overlaps). The features on
  either side of the gene unless those positions exactly
  match the first or last base covered by the gene.
  Therefore, genes falling between two features will at
  least cover two features. This is similar to
  rangeSampleMeans, but it checks the subject positions for
  being sorted and not being NA and also treats them as
  doubles, not ints. Range bounding performed by the
  boundingIndices function.
}
\examples{
data(genoset)
  my.genes = RangedData( ranges=IRanges(start=c(35e6,128e6),end=c(37e6,129e6),names=c("HER2","CMYC")), space=c("chr17","chr8"), universe="hg19")
  rangeSampleMeans( my.genes, genoset.ds, "lrr" )
}
\seealso{
  Other "range summaries": \code{\link{boundingIndices}},
  \code{\link{boundingIndices2}},
  \code{\link{boundingIndicesByChr}},
  \code{\link{rangeColMeans}}
}

