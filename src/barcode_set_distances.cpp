// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "distance.h"
#include "sequencelevenshteindistance.h"
#include "levenshteindistance.h"
#include "hammingdistance.h"

#include "create_distance_func.h"

using namespace Rcpp;

// [[Rcpp::export(".barcode_set_distances")]]
NumericMatrix barcode_set_distances(std::string metric, std::vector< std::string > barcodes, unsigned int cores, const unsigned int cost_sub, const unsigned int cost_indel) {

#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif

  size_t n = barcodes.size();

  NumericMatrix distances(n,n);

  // Convert strings to sequences
  std::vector<Sequence> seqs;
  seqs.reserve(n);

  for (size_t i(0); i < n; i++) {
    seqs.push_back(Sequence(barcodes[i]));
  }

  boost::shared_ptr< Distance > dist = create_distance_func(metric, cost_sub, cost_indel);

  for (size_t i(0); i < n; i++) {
    // Set the diagonal to zero
    distances(i,i) = 0;
    for (size_t j(0); j < i; j++) {
      // Set the symmetrical values of both sequences
      distances(i,j) = dist->distance(seqs[i],seqs[j]);
      distances(j,i) = distances(i,j);
    }
  }

  return distances;
}

