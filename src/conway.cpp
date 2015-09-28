// in conway.cpp

#include "conway.h"
#include "distance.h"
#include "sequence.h"
#include "helpers.h"

#include <vector>
#include <algorithm>

#include <boost/foreach.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

std::vector<Sequence> Conway::close(const std::vector<Sequence> &seed, const std::vector<Sequence> &pool, boost::shared_ptr<Distance> dist, const unsigned int min_dist, const size_t n) {
  std::vector<Sequence> code(seed);

  // walk through every sequence of the pool
  BOOST_FOREACH( const Sequence &seq , pool) {
#ifdef _OPENMP
    if (omp_get_thread_num() == 0) // only in master thread!
#endif
      if (check_interrupt()) {
        // Not in a for loop
        return(std::vector<Sequence>());
      }

    if (dist->is_seq_insertable(code, seq, n, min_dist))
      code.push_back(seq);
  }

  return(code);
}

std::vector<Sequence> Conway::close(const std::vector<Sequence> &seed, const std::vector<Sequence> &pool, boost::shared_ptr<Distance> dist, const unsigned int min_dist, const size_t n, const std::vector<Sequence> &appendix) {
  std::vector<Sequence> code(seed);

  std::vector<Sequence> appended_seed;

  BOOST_FOREACH( const Sequence &seq, seed) {
    BOOST_FOREACH( const Sequence &a,  appendix) {
      appended_seed.push_back(seq.append(a));
    }
  }

  // walk through every sequence of the pool
  BOOST_FOREACH( const Sequence &seq,  pool) {
    bool addable = true;

    // Check for interrupt
#ifdef _OPENMP
    if (omp_get_thread_num() == 0) // only in master thread!
#endif
      if (check_interrupt()) {
        // Not in a for loop
        return(std::vector<Sequence>());
      }
    BOOST_FOREACH( const Sequence &a, appendix) {
      Sequence s_a = seq.append(a);
      if (!dist->is_seq_insertable(appended_seed, s_a, n, min_dist)) {
        addable = false;
        break;
      }
    }
    if (addable) {
      code.push_back(seq);
      BOOST_FOREACH(const Sequence &a, appendix) {
        appended_seed.push_back(seq.append(a));
      }
    }
  }

  return(code);
}

