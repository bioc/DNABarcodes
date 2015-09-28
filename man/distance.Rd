\name{distance}
\alias{distance}
\title{ Calculate distance between two barcodes. }
\description{

The function calculates the distance between two barcodes.
The user may choose one of several distance metrics (\code{"hamming"}, \code{"seqlev"},
\code{"levenshtein"}, \code{"phaseshift"}).

}
\usage{
distance(sequence1, sequence2, metric=c("hamming","seqlev","levenshtein", "phaseshift"), cost_sub=1, cost_indel=1) 
}
\arguments{
	\item{sequence1}{ The first sequence (a string)}
	\item{sequence2}{ The second sequence (a string)}
	\item{metric}{The distance metric which should be calculated.}
	\item{cost_sub}{The cost weight given to a substitution.}
	\item{cost_indel}{The cost weight given to insertions and deletions.}
}
\value{ The distance between the two sequences. }
\examples{
distance("AGGT", "TTCC")
distance("AGGT", "TTCC", metric="seqlev")
}
\references{

Buschmann, T. and Bystrykh, L. V. (2013) Levenshtein error-correcting barcodes for multiplexed DNA sequencing. BMC bioinformatics, 14(1), 272.  Available from \url{http://www.biomedcentral.com/1471-2105/14/272}.

Levenshtein, V. I. (1966). Binary codes capable of correcting deletions, insertions and reversals. In Soviet physics doklady (Vol. 10, p. 707).

Hamming, R. W. (1950). Error detecting and error correcting codes. Bell System technical journal, 29(2), 147-160.

}
