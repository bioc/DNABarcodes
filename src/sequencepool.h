// sequencepool.h

#ifndef __SEQUENCEPOOL_H_INCLUDED__ 
#define __SEQUENCEPOOL_H_INCLUDED__

#include <vector>
#include <memory>

#include "distance.h"

#include "Rcpp.h"

using namespace Rcpp;

class SequencePool { 

  private:
    SequencePool(); // no default constructor

  public:
    static std::vector<Sequence> generate(const size_t, const std::vector< std::string >&, const bool, const bool, const bool);
    static std::vector<Sequence> generate(const size_t, const bool, const bool, const bool);
};

#endif 
