\name{space}
\alias{[}
\alias{elementLengths}
\alias{end}
\alias{loadGC}
\alias{names}
\alias{ranges}
\alias{space}
\alias{start}
\title{Get space factor for GenoSet}
\arguments{
  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{i}{character, RangedData, RangesList, logical,
  integer}

  \item{j}{character, RangedData, RangesList, logical,
  integer}

  \item{k}{character or integer}

  \item{drop}{logical drop levels of space factor?}

  \item{...}{additional subsetting args}

  \item{object}{A GenoSet object or derivative}

  \item{expand}{numeric, expand each feature location by
  this many bases on each side}

  \item{bsgenome,}{sequence db object from BSgenome (e.g.
  Hsapiens)}
}
\value{
  factor

  integer

  integer

  character

  character

  character

  An updated object, with GC percentage information added
  to the locData slot.
}
\description{
  locData slot holds a RangedData, which keeps the
  chromsome of each feature in a factor names 'space'.

  locData slot holds a RangedData.

  locData slot holds a RangedData.

  Get chromosome names

  Get ranges from locData slot

  Get elementLengths from locData slot

  Local GC content can be used to remove GC artifacts from
  copynumber data see Diskin, 2008). GC% column will be
  added to the feature data.  The dataset may be truncated
  to remove positions without GC information.  GC data are
  accessible with locData(). Uses a cool BSgenome trick
  from Michael Lawrence. This takes 5.6 hours for 2Mb
  windows on 2.5M probes, so look for some custom C in
  future releases.
}
\details{
  Get chromosome names, which are the names of the locData
  slot.

  Get ranges from locData slot

  Get elementLengths from locData slot
}
\examples{
data(genoset)
space(genoset.ds)
start(genoset.ds)
end(genoset.ds)
names(genoset.ds)
ranges(genoset.ds) # Returns a RangesList
elementLengths(genoset.ds) # Returns the number of probes per chromosome
data(genoset)
  genoset.ds[1:5,2:3]  # first five probes and samples 2 and 3
  genoset.ds[ , "K"]  # Sample called K
  rd = RangedData(ranges=IRanges(start=seq(from=15e6,by=1e6,length=7),width=1),names=letters[8:14],space=rep("chr17",7))
  genoset.ds[ rd, "K" ]  # sample K and probes overlapping those in rd, which overlap specifed ranges on chr17
}
\author{
  Peter M. Haverty

  Peter M. Haverty

  Peter M. Haverty

  Peter Haverty

  Peter Haverty

  Peter Haverty

  Peter M. Haverty
}

