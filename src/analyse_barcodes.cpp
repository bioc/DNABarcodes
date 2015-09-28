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

#include <boost/shared_ptr.hpp>

using namespace Rcpp;

// [[Rcpp::export(".analyse_barcodes")]]
std::vector<double> analyse_barcodes(const std::string metric, const std::vector< std::string > barcodes, const unsigned int cores, const unsigned int cost_sub, const unsigned int cost_indel ) {

#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif

  size_t n = barcodes.size();

  size_t n_pairs = (n * n - n) / 2; // will always be without remainder

  // Convert strings to sequences
  std::vector<Sequence> seqs;
  for (size_t i(0); i < n; i++) {
    seqs.push_back(Sequence(barcodes[i]));
  }

  boost::shared_ptr< Distance > dist = create_distance_func(metric, cost_sub, cost_indel);

  std::vector<double> return_distances(n_pairs);

  size_t pair_iterator = 0;

  for (size_t i(0); i < n; i++) {
    for (size_t j(0); j < i; j++) {
      return_distances[pair_iterator] = dist->distance(seqs[i],seqs[j]);
      pair_iterator++;
    }
  }

  return return_distances;
}

