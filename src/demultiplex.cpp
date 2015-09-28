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
#include "phaseshiftdist.h"

#include "create_distance_func.h"

using namespace Rcpp;

// [[Rcpp::export(".demultiplex")]]
Rcpp::DataFrame demultiplex(std::vector< std::string > barcodes, std::vector< std::string > reads, std::string metric, const unsigned int cost_sub, const unsigned int cost_indel) {
  // Get the right distance algorithm
  boost::shared_ptr< Distance > dist = create_distance_func(metric, cost_sub, cost_indel);

  return dist->demultiplex(barcodes,reads);
}
