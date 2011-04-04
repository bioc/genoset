\name{universe}
\alias{universe,-method}
\alias{universe}
\title{Get and set the genome universe annotation.}
\description{Genome universe for locData}
\details{\code{universe,-method}: The genome positions of the features in locData. The UCSC notation (e.g. hg18, hg19, etc.) should be used.

\code{universe<-,-method}: Set genome universe


}
\value{\code{universe,-method}: character, e.g. hg19

\code{universe<-,-method}: A GenoSet object


}
\author{Peter M. Haverty}
\arguments{\item{x}{GenoSet}
\item{value}{character, new universe string, e.g. hg19}
}
\examples{data(genoset)
universe(genoset.ds)
universe(genoset.ds) = "hg19"}
\alias{universe<-,-method}
\alias{universe<-}

