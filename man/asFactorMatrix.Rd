\name{asFactorMatrix}
\alias{asFactorMatrix}
\title{Make factor matrix from character matrix...}
\usage{asFactorMatrix(object, levels)
}
\description{Make factor matrix from character matrix}
\details{Make factor matrix from character matrix for use with convertToBigMatrix.
Makes an integer matrix with levels since as.big.matrix would make a
factor matrix into a 1D object for some reason. Character matrices should
be converted to factors with explicit levels as huge matrices are likely
too big to unique.

Caution: use asFactorMatrix on matrices already in an eSet.  The eSet constructor will
apparently wipe out the levels.}
\value{factor with dimensions matching object}
\author{Peter M. Haverty \email{phaverty@gene.com}}
\arguments{\item{object}{matrix of characters}
\item{levels}{character}
}

