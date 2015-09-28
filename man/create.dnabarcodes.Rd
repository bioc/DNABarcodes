\name{create.dnabarcodes}
\alias{create.dnabarcodes}
\alias{create.dnabarcodes.clique}
\alias{create.dnabarcodes.conway}
\alias{create.dnabarcodes.sampling}
\alias{create.dnabarcodes.ashlock}
\title{ Create a set of DNA barcodes using one of several heuristic methods.}
\description{

Creates a DNA barcode set that has error correction and detection properties.
The function uses one of four different heuristics to generate the set (clique,
conway, sampling, and ashlock) and one of three different distance metrics between
individual barcodes (hamming, seqlev, and levenshtein).

For heuristics, inexperienced users should try "conway" and then "ashlock".

The Hamming Distance (\code{metric="hamming"}) allows the correction/detection of
substitutions. The Sequence Levenshtein distance (\code{metric="seqlev"}) allows the
correction/detection of insertions, deletions, and substitutions in barcodes in
DNA context. The Levenshtein distance should not be used except by experienced
users (see Details).

The functions \code{create.dnabarcodes.conway},
\code{create.dnabarcodes.clique}, \code{create.dnabarcodes.sampling}, and
\code{create.dnabarcodes.ashlock} provide a shortcut to generating barcode sets
based on one of the available heuristics and use better default parameters for
some.
}
\usage{
create.dnabarcodes(n, dist=3, metric=c("hamming","seqlev","levenshtein", "phaseshift"), 
                   heuristic=c("conway", "clique", "sampling", "ashlock"), 
                   filter.triplets=TRUE, filter.gc=TRUE, filter.self_complementary=TRUE, 
                   pool = character(), iterations=100, population=200, 
                   cores=detectCores()/2, use_cache = FALSE,
                   cost_sub = 1, cost_indel = 1) 

create.dnabarcodes.conway(n, heuristic="conway", \dots)

create.dnabarcodes.clique(n, heuristic="clique", \dots)

create.dnabarcodes.sampling(n, heuristic="sampling", iterations=20000, \dots)

create.dnabarcodes.ashlock(n, heuristic="ashlock", iterations=100, population=200, \dots)
}
\details{

Different heuristics produce different results in different time. New users
should first try "conway" and then "ashlock".  

The heuristics "conway" and "clique" produce fast results but are not nearly as
good as the heuristics "sampling" and "ashlock". The clique heuristic is a bit
slower and needs more memory than the Conway heuristic because it first
constructs a graph representation of the pool.

The heuristic "ashlock" is assumed to produce the best heuristic results after
a reasonable number of iterations with a good population size.  

Distance metrics are the mathematical fairy dust that make the error correction
and detection of the barcode sets possible. Different metrics and different
distances allow different error corrections/detections.

A high enough Hamming Distance (\code{metric = "hamming"}) allows the
correction/detection of substitutions. Due to the ignorance of insertions and
deletions, any changes to the length of the barcode as well as DNA context are
ignored, which makes the Hamming distance a simple choice.

A high enough Sequence Levenshtein distance (\code{metric = "seqlev"}) allows the
correction/detection of insertions, deletions, and substitutions in scenarios
where the barcode was attached to a DNA sequence. Your sequence read, coming
from a Illumina, Roche, PacBio NGS machine of your choice, should then start
with that barcode, followed by the insert or an adapter or some random base
calls.

A high enough Levenshtein distance (\code{metric = "levenshtein"}) allows the
correction/detection of insertions, deletions, and substitutions in scenarios
where the barcode was not attached anywhere, respective where the exact outline
of the barcode is known. This is as far as we know in no current NGS technology
the case. Do not use this distance metric, except you know what you are doing.

The number of error corrections/detections of the code depends of the enforced
distance \code{dist}. If all conditions are correct, a barcode set with an
enforced distance \code{dist} can correct \code{k} errors, if 

\deqn{ k \le floor(\frac{dist - 1}{2})}{k <= floor((dist-1)/2)}

The detection of \code{k} errors is possible, if

\deqn{ k \le dist - 1}{k <= dist - 1}

The advantage of this function over already available barcode sets in the
scientific community is the ability to flexibly generate new barcode sets of
different properties. For example, the function can use a pre-existing barcode
library as a candidate set for a better barcode set. In another example, a
higher distance (e.g., \code{dist = 4}) is used. Such a parameter setting would
possibly increase the error detection property of the code as well as the
average barcode distance, increasing the probability of guessing a barcode
during demultiplexing.

The heuristics for the generation of barcodes are as follows:

The Conway heuristic (\code{heuristic = "conway"}, named after John Conway)
starts with an empty set of barcodes, goes through the list of candidate
barcodes (the \code{pool}) in lexicographical order and adds each candidate
barcode to the initial set if the distance if the candidate barcode to each
barcode in the intial set is at least \code{d >= dist}.

The Clique heuristic (\code{heuristic = "conway"}) first generates a graph
representation of the pool. Each barcode in the pool is a node of the graph and
two barcodes/nodes are connected undirectionally if their distance is at least
\code{d >= dist}. The barcode set problem is now reduced to finding the
maximal clique in this graph. Because that problem is also computationally
infeasible, we use the heuristic clique algorithm of Pattabiraman et al.

The sampling heuristic (\code{heuristic = "sampling"}) extends the principle of
the Conway heuristic. Instead of starting with an empty initial set, we
generate small random sets of barcodes as initial sets (the so called seeds).
Those seeds are then "closed" using the Conway method. The size of the seed is
fixed to three barcodes. The number of random seeds is given by the parameter
\code{iterations}, hence that often a Conway closure is calculated.

Finally, the Ashlock heuristic (\code{heuristic = "ashlock"}, named after
Daniel Ashlock) extends the sampling heuristic by adding an evolutionary
algorithm. A population of random seeds is generated only for the first
iteration. Each seed is then closed using the Conway method. The size of the
barcode set after closure defines the \emph{fitness} of that seed. For the next
iteration, succesful seeds (with a higher fitness) are
cloned and slightly mutated (some barcodes in the seed are replaced with a
random new barcode). Those changed seeds are now closed again and their
respective fitness calculated. In the first iteration, as many instances of the
Conway closure are calculated as there are seeds. In the second round, only
half of the seeds (the changed ones) are calculated. Therefore, the total
number of calculated Conway closures is \code{iterations + population/2 *
(iterations - 1)}.  } \arguments{ \item{n}{The length of the barcodes (should
be smaller than 14). No default.}

	\item{dist}{The minimal distance between barcodes that shall be kept. Default is 3.}

	\item{metric}{The distance metric that is used to calculate and keep
          the distance between barcodes. Default is "hamming"}

	\item{heuristic}{The heuristic algorithm to generate the barcode set.
          Available are "conway", "clique", "sampling", and "ashlock". The default is
          "conway". }

	\item{filter.triplets}{Should sequences that contain at least three
          repeated equal bases (e.g., AAA, TTT, CCC, or GGG) be filtered out?}

	\item{filter.gc}{Should sequences that have an unbalanced ratio of
          bases G or C versus A or T be filtered out?}

	\item{filter.self_complementary}{Should self complementary sequences be
          filtered out?}

	\item{pool}{An optional set of candidate sequences for the DNA barcode
sets. If no pool is given, the maximum possible pool is generated internally in
the function according to specifications (length \code{n}, filtered based on
\code{filter.triplets}, \code{filter.gc}, \code{filter.self_complementary}).
If the pool is given, only sequences from the pool are used as possible
barcodes (but still filtered according to filtering parameters of this
function).  It is save to leave this option unset.}

	\item{iterations}{In case of \code{heuristic = "sampling"}: The number of
samples that are tested for maximal size. In case of \code{heuristic = "ashlock"}: The
number of iterations of the genetic algorithm that are conducted. Not used in any
other case.}

	\item{population}{Only used for \code{heuristic = "ashlock"}: The
number of chromosomes of the genetic algorithm that are tested. Note: For
\code{heuristic = "ashlock"}, the number of barcode sets that are tested is
population + population/2 * (iterations-1).}

	\item{cores}{The number of cores (CPUs) that will be used for parallel
(openMP) calculations.}
	
        \item{use_cache}{Shall the distances between each candidate barcode of the pool be calculated in advance? In many cases this increases speed but needs a lot of memory. When in doubt, set to FALSE.}
	
        \item{cost_sub}{The cost weight given to a substitution.}

	\item{cost_indel}{The cost weight given to insertions and deletions.}

	\item{\dots}{Arguments passed on to \code{create.dnabarcodes}.}

}
\value{ A vector of characters, representing the DNA barcode set. }
\examples{
# Create barcodes of length 5 and minimal Hamming distance of 3:
# (these barcodes can correct up to 1 substitution mutation)
create.dnabarcodes(5) # 30 barcodes
create.dnabarcodes.ashlock(5) # Up to 48 barcodes

# Create barcodes of length 5 with minimal SeqLev distance of 3:
# (barcodes can correct up to 1 insertion/deletion/substitution
#  in DNA context)
create.dnabarcodes(5, metric="seqlev") # 8 barcodes
create.dnabarcodes.ashlock(5, metric="seqlev") # Up to 13

# Create Seq-Lev-barcodes without filtering
length(create.dnabarcodes.ashlock(5, metric="seqlev", filter.triplets=FALSE, filter.gc=FALSE, filter.self_complementary=FALSE)) # Up to 15

# Create pool, apply additional (useless) filters, create barcode set
pool <- create.pool(5)
length(pool) # 592
pool <- pool[grep("AT",pool,invert=TRUE)]
length(pool) # 468
create.dnabarcodes(5,pool=pool)
}
\seealso{
\code{\link{create.pool}} to get the pool of barcode candidates that are used in this function internally.
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

