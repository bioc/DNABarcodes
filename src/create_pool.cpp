// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#include <boost/foreach.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "sequence.h"
#include "sequencepool.h"
#include "helpers.h"

using namespace Rcpp;

// [[Rcpp::export(".create_pool")]]
std::vector< std::string > create_pool(unsigned long int n, bool filter_triplets, bool filter_gc, bool filter_self_complementary, unsigned int cores) {

#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif

  std::vector<Sequence> pool = SequencePool::generate(n, filter_triplets, filter_gc, filter_self_complementary);

  std::vector< std::string > str_pool;
  str_pool.reserve(pool.size());

  BOOST_FOREACH(Sequence seq, pool) {
#ifdef _OPENMP
          if (omp_get_thread_num() == 0) // only in master thread!
#endif
            if (check_interrupt()) {
              // Not in a for loop
              return(std::vector< std::string >());
            }
          str_pool.push_back(seq.asString());
  }

  return str_pool;
}

