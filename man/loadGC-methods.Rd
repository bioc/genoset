\docType{methods}
\name{loadGC}
\alias{loadGC}
\alias{loadGC,GenoSet,numeric,BSgenome-method}
\alias{loadGC,RangedData,numeric,BSgenome-method}
\title{Load local GC percentage around features}
\arguments{
  \item{object}{A GenoSet object or derivative}

  \item{expand}{numeric, expand each feature location by
  this many bases on each side}

  \item{bsgenome,}{sequence db object from BSgenome (e.g.
  Hsapiens)}
}
\value{
  An updated object, with GC percentage information added
  to the locData slot.
}
\description{
  Local GC content can be used to remove GC artifacts from
  copynumber data see Diskin, 2008). GC% column will be
  added to the feature data.  The dataset may be truncated
  to remove positions without GC information.  GC data are
  accessible with locData(). Uses a cool BSgenome trick
  from Michael Lawrence. This takes 5.6 hours for 2Mb
  windows on 2.5M probes, so look for some custom C in
  future releases.
}
\author{
  Peter M. Haverty
}
\seealso{
  Other "gc content": \code{\link{gcCorrect}}
}

