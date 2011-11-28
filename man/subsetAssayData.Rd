\name{subsetAssayData}
\alias{subsetAssayData}
\title{Subset assayData}
\usage{subsetAssayData(orig, i, j, ..., drop=FALSE)}
\description{Subset or re-order assayData}
\details{Subset or re-order assayData locked environment, environment, or list. Shamelessly stolen
from "[" method in Biobase version 2.8 along with guts of assayDataStorageMode()}
\value{assayData data structure}
\author{Peter M. Haverty}
\arguments{\item{orig}{assayData environment}
\item{i}{row indices}
\item{j}{col indices}
\item{...}{Additional args to give to subset operator}
\item{drop}{logical, drop dimensions when subsetting with single value?}}
\examples{data(genoset)
ad = assayData(genoset.ds)
small.ad = subsetAssayData(ad,1:5,2:3)}
