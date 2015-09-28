\name{barcode.set.distances}
\alias{barcode.set.distances}
\title{ Calculate distances between each barcode pair of a barcode set. }
\description{

The function calculates the distance between each pair of a set of barcodes.
The user may choose one of several distance metrics (\code{"hamming"}, \code{"seqlev"},
\code{"levenshtein"})).

}
\usage{
barcode.set.distances(barcodes, metric=c("hamming", "seqlev", "levenshtein"), cores=detectCores()/2, cost_sub = 1, cost_indel = 1)	
}
\arguments{
	\item{barcodes}{A set of barcodes (as vector of characters)}
	\item{metric}{The distance metric which should be calculated.}
	\item{cores}{The number of cores (CPUs) that will be used for parallel (openMP) calculations.}
	\item{cost_sub}{The cost weight given to a substitution.}
	\item{cost_indel}{The cost weight given to insertions and deletions.}
}
\value{ A symmetric \code{Matrix} of distances between each pair of barcodes with zeros on the main diagonal. }
\examples{
barcodes <- c("AGGT", "TTCC", "CTGA", "GCAA")
barcode.set.distances(barcodes)
barcode.set.distances(barcodes,metric="seqlev")
}
\details{
The primary purpose of this function is the analysis of barcode sets. Seeing
the individual paired barcode distances helps to understand which pairings are
exceptionally similar and which barcodes have a smaller or higher average
distance to other barcodes.

Details if the distance metrics can be found in the man page of \code{\link{create.dnabarcodes}}.
}
\seealso{
\code{\link{analyse.barcodes}}
}
\references{

Buschmann, T. and Bystrykh, L. V. (2013) Levenshtein error-correcting barcodes for multiplexed DNA sequencing. BMC bioinformatics, 14(1), 272.  Available from \url{http://www.biomedcentral.com/1471-2105/14/272}.

Levenshtein, V. I. (1966). Binary codes capable of correcting deletions, insertions and reversals. In Soviet physics doklady (Vol. 10, p. 707).

Hamming, R. W. (1950). Error detecting and error correcting codes. Bell System technical journal, 29(2), 147-160.

}
