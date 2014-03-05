\name{segs2RangedData}
\alias{segs2RangedData}
\title{Make a RangedData from segments}
\usage{
segs2RangedData(segs)
}
\arguments{
  \item{segs}{data.frame, like from segment in DNAcopy or
  segTable}
}
\value{
RangedData
}
\description{
Starting from a data.frame of segments, like from CBS and
segTable, organize as a RangedData. Label data "score", so
it can easily be made into various genome browser formats
using rtracklayer.
}
\author{
Peter M. Haverty \email{phaverty@gene.com}
}
\seealso{
Other "segmented data": \code{\link{bounds2Rle}};
\code{\link{rangeSegMeanLength}},
\code{\link{rangeSegMeanLength,GRanges,data.frame-method}},
\code{\link{rangeSegMeanLength,GRanges,data.frame-method}},
\code{\link{rangeSegMeanLength,GRanges,list-method}},
\code{\link{rangeSegMeanLength,GRanges,list-method}};
\code{\link{runCBS}}, \code{\link{runCBS}};
\code{\link{segPairTable}},
\code{\link{segPairTable,DataFrame,DataFrame-method}},
\code{\link{segPairTable,DataFrame,DataFrame-method}},
\code{\link{segPairTable,Rle,Rle-method}},
\code{\link{segPairTable,Rle,Rle-method}};
\code{\link{segs2Granges}}; \code{\link{segs2Rle}};
\code{\link{segs2RleDataFrame}}; \code{\link{segTable}},
\code{\link{segTable,DataFrame-method}},
\code{\link{segTable,DataFrame-method}},
\code{\link{segTable,Rle-method}},
\code{\link{segTable,Rle-method}}

Other segments: \code{\link{segs2RleDataFrame}}
}

