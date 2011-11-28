\name{boundingIndices2}
\alias{boundingIndices2}
\title{Find indices of features bounding a set of chromsome ranges/genes...}
\usage{boundingIndices2(starts, stops, positions, offset)}
\description{Find indices of features bounding a set of chromsome ranges/genes}
\details{This function is similar to findOverlaps but it guarantees at least two features will be
covered. This is useful in the case of finding features corresponding to a set of genes.
Some genes will fall entirely between two features and thus would not return any ranges
with findOverlaps. Specifically, this function will find the indices of the features
(first and last) bounding the ends of a range/gene (start and stop) such that
first <= start <= stop <= last. Equality is necessary so that multiple conversions between
indices and genomic positions will not expand with each conversion. This function uses
findIntervals, which is for k queries and n features is O(k * log(n)) generally and
~O(k) for sorted queries. Therefore will be dramatically faster for sets of query genes
that are sorted by start position within each chromosome.  This should give performance
for k genes and n features that is ~O(k) for starts and O(k * log(n)) for stops and
~O(k * log(n)) overall.  Ranges/genes that are outside the range of feature positions will
be given the indices of the corresponding first or last index rather than 0 or n + 1 so
that genes can always be connected to some data.}
\value{integer matrix of 2 columms for start and stop index of range in data}
\author{Peter M. Haverty}
\arguments{\item{starts}{numeric or integer, first base position of each query range}
\item{stops}{numeric or integer, last base position of each query range}
\item{positions}{Base positions in which to search}
\item{offset}{integer, value to add to all returned indices. For the case where positions represents a portion of some larger array (e.g. a chr in a genome)}}
\examples{starts = seq(10,100,10)
boundingIndices2( starts=starts, stops=starts+5, positions = 1:100 )}
