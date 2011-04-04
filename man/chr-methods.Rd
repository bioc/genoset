\name{chr-methods}
\alias{chr-methods}
\alias{chr}
\title{Look up chromosome for each feature}
\description{Chromsome name for each feature}
\details{\code{chr-methods}: Get chromosome name for each feature.  Returns character, not the factor 'space'.

}
\value{\code{chr-methods}: character vector of chromosome positions for each feature

}
\author{Peter Haverty}
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
chr(gs)  # c("chr1","chr1","chr1","chr1","chr3","chr3","chrX","chrX","chrX","chrX")
chr(locData(gs))  # The same}
\alias{chr,-method}
\alias{chr}
\alias{chr,-method}
\alias{chr}

