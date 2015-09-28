\name{analyse.barcodes}
\alias{analyse.barcodes}
\title{Analyse Sets of Barcode}
\description{
The function analyses the properties of sets of (potential) barcodes. For
various metrics, minimal, maximal and average distances are calculated and
hints on error correction capabilities of the code are given.
}
\usage{
analyse.barcodes(barcodes, metric = c("hamming", "seqlev", "levenshtein"), cores=detectCores()/2, cost_sub = 1, cost_indel = 1)
}
\arguments{

	\item{barcodes}{A vector of characters that represent the barcodes. All barcodes must be of equal length and consist only of letters A, C, G, and T. Lower case letters are allowed but do not make a difference.}

	\item{metric}{A vector of one or more metric names whose distances shall be calculated for the barcode set. Default is to use all of them.}
	
	\item{cores}{The number of cores (CPUs) that will be used for parallel (openMP) calculations.}

	\item{cost_sub}{The cost weight given to a substitution.}

	\item{cost_indel}{The cost weight given to insertions and deletions.}

}
\value{ A data frame of properties of the barcode set with the following meanings:

Columns: The first column contains a description of the barcode set properties. Each next column names the metric.

Rows: "Mean Distance", "Median Distance", "Minimum Distance", "Maximum Distance", "Guaranteed Error Correction", "Guaranteed Error Detection" 

}
\examples{
barcodes <- c("ACG", "CGT", "TGC")
analyse.barcodes(barcodes)

##
##                 Description  hamming   seqlev levenshtein
##               Mean Distance 2.666667 1.666667    2.333333
##             Median Distance 3.000000 2.000000    2.000000
##            Minimum Distance 2.000000 1.000000    2.000000
##            Maximum Distance 3.000000 2.000000    3.000000
## Guaranteed Error Correction 0.000000 0.000000    0.000000
##  Guaranteed Error Detection 1.000000 0.000000    1.000000

}
\seealso{
\code{\link{barcode.set.distances}}
}
\references{

Buschmann, T. and Bystrykh, L. V. (2013) Levenshtein error-correcting barcodes for multiplexed DNA sequencing. BMC bioinformatics, 14(1), 272.  Available from \url{http://www.biomedcentral.com/1471-2105/14/272}.

Levenshtein, V. I. (1966). Binary codes capable of correcting deletions, insertions and reversals. In Soviet physics doklady (Vol. 10, p. 707).

Hamming, R. W. (1950). Error detecting and error correcting codes. Bell System technical journal, 29(2), 147-160.

}
