\name{DNABarcodes-package}
\alias{DNABarcodes-package}
\alias{DNABarcodes}
\docType{package}
\title{
Create and analyse DNA barcode sets that are capable of error correction.
}
\description{
The package offers a function to create DNA barcode sets capable of correcting
substitution errors or insertion, deletion, and substitution errors. Existing
barcodes can be analysed regarding their minimal, maximal and average distances
between barcodes. Finally, reads that start with a (possibly mutated) barcode
can be demultiplexed, i.e. assigned to their original reference barcode.
}
\details{
\tabular{ll}{
Package: \tab DNABarcodes\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-07-23\cr
License: \tab GPL-2\cr
}

The function \code{\link{create.dnabarcodes}} creates a set of barcodes of
equal length that satisfies some wished criteria regarding error correction.

After sequencing the DNA/RNA material, the researcher will have a set of reads
that start with a (possibly mutated) barcode. For Illumina HiSeq, this is the
index read. For PacBio, this is the read itself (with some other
complications). The function \code{\link{demultiplex}} can then be used to
assign reads to their original reference barcodes. \code{\link{demultiplex}}
will correct mutations in a best-effort way.

Existing sets of barcodes (e.g. supplied by a manufacturer) can be analysed
with functions \code{\link{analyse.barcodes}} and
\code{\link{barcode.set.distances}}.

The advantage of this package over using already available barcode sets in the
scientific community is the ability to flexibly generate new barcode sets of
different properties. For example, \code{\link{create.dnabarcodes}} can use a
pre-existing barcode library as a candidate set for a better barcode set. In
another example, a higher distance (e.g., \code{dist = 4}) can be used. Such a
parameter setting would possibly increase the error detection property of the
code as well as the average barcode distance, increasing the probability of
guessing a barcode during demultiplexing.

}
\author{
Tilo Buschmann (tilo.buschmann.ac@gmail.com)
}
\references{

Buschmann, T. and Bystrykh, L. V. (2013) Levenshtein error-correcting barcodes for multiplexed DNA sequencing. BMC bioinformatics, 14(1), 272.  Available from \url{http://www.biomedcentral.com/1471-2105/14/272}.

Levenshtein, V. I. (1966). Binary codes capable of correcting deletions, insertions and reversals. In Soviet physics doklady (Vol. 10, p. 707).

Hamming, R. W. (1950). Error detecting and error correcting codes. Bell System technical journal, 29(2), 147-160.

Conway, J. and Sloane, N. (1986) Lexicographic codes: error-correcting codes from game theory. Information Theory, IEEE Transactions on, 32(3), 337-348.

Pattabiraman, B., Patwary, M. M. A., Gebremedhin, A. H., Liao, W. K. and Choudhary, A. (2013) Fast algorithms for the maximum clique problem on massive sparse graphs. In Algorithms and Models for the Web Graph (pp. 156-169).  Springer International Publishing.  

Ashlock, D., Guo, L. and Qiu, F. (2002) Greedy closure evolutionary algorithms.  In Computational Intelligence, Proceedings of the World on Congress on (Vol. 2, pp. 1296-1301). IEEE.

Brouwer, A. E., Shearer, L. B. and Sloane, N. I. A. (1990) A new table of constant weight codes. In IEEE Trans Inform Theory.

}
\concept{Error Correction}
\concept{DNA Barcodes}
\concept{Sample Tags}
\concept{Levenshtein Distance}
\concept{Hamming Distance}
\concept{Next Generation Sequencing}
\concept{Multiplexing}
\concept{Demultiplexing}

\examples{
# Create Sequence Levenshtein Barcodes with the default heuristic
dnabarcodes1 <- create.dnabarcodes(5, metric="seqlev")

# Create Sequence Levenshtein Barcodes with a better, but slower heuristic
dnabarcodes2 <- create.dnabarcodes(5, metric="seqlev", heuristic="ashlock")
}
