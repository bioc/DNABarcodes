\name{demultiplex}
\alias{demultiplex}
\title{ Demultiplex a set of reads. }
\description{

The function \code{demultiplex} takes a set of reads that start with a barcode
and assigns those reads to a reference barcode while possibly correcting
errors.

The correct metric should be used, with \code{metric = "hamming" } to correct
substitution errors and \code{metric = "seqlev" } to correct insertion,
deletion, and substitution errors.
}
\note{
Do not try to correct errors in barcodes that were not systematically constructed for such a correction. To create such a barcode set, have a look into function \code{\link{create.dnabarcodes}}.
}
\usage{
demultiplex(reads, barcodes, metric=c("hamming","seqlev","levenshtein","phaseshift"), cost_sub = 1, cost_indel = 1)
}
\arguments{

	\item{reads}{The reads coming from your sequencing machines that start
with a barcode. For \code{metric = "seqlev"} please provide some context after
the (supposed) barcode, at least as many bases as errors that you want to
correct.}

	\item{barcodes}{The reference barcodes that you used during library
preparation and that you want to correct in your reads.}

	\item{metric}{The distance metric to be used to assign reads to
reference barcodes.}
	
        \item{cost_sub}{The cost weight given to a substitution.}

	\item{cost_indel}{The cost weight given to insertions and deletions.}

}
\details{

Reads are matched to their correct reference barcodes by calculating the
distances between each read and each reference barcode. The reference barcode
with the smallest distance to the read is assumed to be the correct original
barcode of that read.

For \code{metric = "hamming"}, only the first \code{n} (with \code{n} being the
length of the reference barcodes) bases of the read are used for these
comparisons and no bases afterwards. Reads with fewer than \code{n} bases
cannot be matched.

For \code{metric = "seqlev"}, the whole read is compared with the reference
barcodes. The Sequence Levenshtein distance was especially developed for
barcodes in DNA context and can cope with ambiguities that stem from changes
to the length of the barcode.

The Levenshtein distance (\code{metric = "levenshtein"}) is largely undefined
in DNA context and should be avoided. The Levenshtein distance only works if
the length both of the reference barcode and the barcode in the read is known.
With possible insertions and deletions, this becomes an unknown. For this
reason, we always calculate the Levenshtein distance between the whole read and
the whole reference barcode without coping with potential side effects.

}

\value{ A vector of reference barcodes of the same length as the input reads.
Each reference barcode is the corrected version of the input barcode.}

\examples{
# Define some barcodes and inserts
barcodes <- c("AGGT", "TTCC", "CTGA", "GCAA")
insert <- 'ACGCAGGTTGCATATTTTAGGAAGTGAGGAGGAGGCACGGGCTCGAGCTGCGGCTGGGTCTGGGGCGCGG'

# Choose and mutate a couple of thousand barcodes
used_barcodes <- sample(barcodes,10000,replace=TRUE)
mutated_barcodes <- unlist(lapply(strsplit(used_barcodes,""), function(x) { pos <- sample(1:length(x),1); x[pos] <- sample(c("C","G","A","T"),1); return(paste(x,collapse='')) } ))

show(setequal(mutated_barcodes, used_barcodes)) # FALSE

# Construct reads (= barcodes + insert)
reads <- paste(mutated_barcodes, insert, sep='')

# Demultiplex
demultiplexed <- demultiplex(reads,barcodes,metric="hamming")

# Show correctness
show(setequal(demultiplexed, used_barcodes)) # TRUE

}
\seealso{
\code{\link{create.dnabarcodes}}, \code{\link{analyse.barcodes}}
}
