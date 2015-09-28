// in sequencepool.cpp

// [[Rcpp::depends(BH)]]
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]

#include <cmath>
#include <algorithm>
#include <memory>

#include "distance.h"
#include "sequence.h"
#include "sequencepool.h"
#include "helpers.h"

#include <boost/foreach.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

std::vector<Sequence> SequencePool::generate(const size_t n, const std::vector< std::string > &str_pool, const bool filter_triplets, const bool filter_gc, const bool filter_self_complementary) {
  std::vector<Sequence> seqs;

  BOOST_FOREACH(std::string str_seq, str_pool) {
    Sequence seq(str_seq);
#ifdef _OPENMP
    if (omp_get_thread_num() == 0) // only in master thread!
#endif
      if (check_interrupt()) {
        // Not in a for loop
        return(std::vector<Sequence>());
      }

    if (!(filter_triplets && seq.containsTriples()) && !(filter_self_complementary && seq.isSelfComplementary()) && !(filter_gc && !seq.isGCContentRight())) {
      seqs.push_back(seq);
    }

  }

  //seqs.shrink_to_fit();
  std::sort(seqs.begin(),seqs.end());

  return(seqs);
}

std::vector<Sequence> SequencePool::generate(const size_t n, const bool filter_triplets, const bool filter_gc, const bool filter_self_complementary) {
  std::vector<Sequence> seqs;

  uint64_t max_val = ((uint64_t) 1 << (n + n));

  for (uint64_t value = 0; value < max_val; value++) {
#ifdef _OPENMP
    if (omp_get_thread_num() == 0) // only in master thread!
#endif
      if (check_interrupt()) {
        // Not in a for loop
        return(std::vector<Sequence>());
      }

    // Go through all bases and convert to new three_bit scheme
    uint64_t new_value = 0;//three-bit value
    for (size_t i = 0; i < n; i++) {
      uint64_t old_base = (value >> (2 * i)) & 0x3;
      uint64_t new_base = Sequence::REAL_BASES[old_base];
      new_value = (new_value << SEQUENCE_BITS_PER_BASE) | new_base;
    }

    Sequence seq(new_value, n);

    // I need to sort this statement out

    if (!(filter_triplets && seq.containsTriples()) && !(filter_self_complementary && seq.isSelfComplementary()) && !(filter_gc && !seq.isGCContentRight())) {
      seqs.push_back(seq);
    }
  }

  //seqs.shrink_to_fit();
  std::sort(seqs.begin(),seqs.end());

  return(seqs);
}

