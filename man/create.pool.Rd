\name{create.pool}
\alias{create.pool}
\title{ Create a pool of barcode candidates.}
\description{

This function creates a vector of sequences of equal length with some filtering
applied. These sequences function as the pool of candidate sequences from which
an actual DNA barcode set may be constructed, for example by passing the pool
as a parameter to the function \code{\link{create.dnabarcodes}}.

Sequences in the pool are constructed as all possible concatenations of n bases
C,G,A, and T with some of them filtered out for various technical reasons.

Users of this package usually do not need to use this function as the function
\code{\link{create.dnabarcodes}} already creates such a pool internally.
However, there are some exceptions: A) The user might want to know the set of
barcode candidates from which the final DNA barcode set was generated. B) The
user might want apply some additional filtering to the pool before constructing
a DNA barcode set.

}
\usage{
create.pool(n, filter.triplets=TRUE, filter.gc=TRUE, filter.self_complementary=TRUE, cores=detectCores()/2)	
}
\arguments{

	\item{n}{The length of the sequences in the pool (should be smaller than 20)}

	\item{filter.triplets}{Should sequences that contain at least three repeated equal bases (e.g., AAA, TTT, CCC, or GGG) be filtered out?}

	\item{filter.gc}{Should sequences that have an unbalanced ratio of bases G or C versus A or T be filtered out?}

	\item{filter.self_complementary}{Should self complementary sequences be filtered out?}
	
	\item{cores}{The number of cores (CPUs) that will be used for parallel (openMP) calculations.}

}
\value{ A vector of characters representing the pool of barcode candidates. }
\examples{
create.pool(4)
length(create.pool(5))
length(create.pool(5, filter.triplets=FALSE, filter.gc=FALSE, filter.self_complementary=FALSE))
#
#
# Create pool, apply additional (useless) filters, create barcode set
pool <- create.pool(5)
length(pool) # 592
pool <- pool[grep("AT",pool,invert=TRUE)]
length(pool) # 468
create.dnabarcodes(5,pool=pool)
}
\seealso{
\code{\link{create.dnabarcodes}}
}

