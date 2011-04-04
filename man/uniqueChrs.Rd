\name{uniqueChrs}
\alias{uniqueChrs-methods}
\alias{uniqueChrs}
\title{Get list of unique chromosome names...}
\description{Get list of unique chromosome names}
\details{\code{uniqueChrs-methods}: Get list of unique chromosome names. A synonym for names().

}
\value{\code{uniqueChrs-methods}: character vector with names of chromosomes

}
\author{Peter M. Haverty}
\arguments{\item{object}{RangedData or GenoSet}
}
\examples{test.sample.names = LETTERS[11:13]
probe.names = letters[1:10]
gs = GenoSet(
locData=RangedData(ranges=IRanges(start=1:10,width=1,names=probe.names),space=c(rep("chr1",4),rep("chr3",2),rep("chrX",4)),universe="hg18"),
cn=matrix(31:60,nrow=10,ncol=3,dimnames=list(probe.names,test.sample.names)),
pData=data.frame(matrix(LETTERS[1:15],nrow=3,ncol=5,dimnames=list(test.sample.names,letters[1:5]))),
annotation="SNP6"
)
uniqueChrs(gs) # c("chr1","chr3","chrX")
uniqueChrs(locData(gs))  # The same}
\alias{uniqueChrs,-method}
\alias{uniqueChrs}

