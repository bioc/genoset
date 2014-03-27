\docType{methods}
\name{start,GenoSet-method}
\alias{dim,GenoSet-method}
\alias{elementLengths,GRanges-method}
\alias{elementLengths,GenoSet-method}
\alias{end,GenoSet-method}
\alias{names,GenoSet-method}
\alias{nrow,GRanges-method}
\alias{start,GenoSet-method}
\alias{width,GenoSet-method}
\title{Get start of location for each feature}
\usage{
\S4method{start}{GenoSet}(x)

\S4method{end}{GenoSet}(x)

\S4method{width}{GenoSet}(x)

\S4method{names}{GenoSet}(x)

\S4method{elementLengths}{GenoSet}(x)

\S4method{elementLengths}{GRanges}(x)

\S4method{nrow}{GRanges}(x)

\S4method{dim}{GenoSet}(x)
}
\arguments{
  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}

  \item{x}{GenoSet}
}
\value{
integer

integer

integer

character

character
}
\description{
Get start of location for each feature

Get end of location for each feature

Get width of location for each feature

Get names of data matrices. For the time being, this is
\code{assayDataElementNames}. This function used to do
\code{chrNames}.

Get elementLengths from locData slot
}
\details{
Get elementLengths from locData slot
}

