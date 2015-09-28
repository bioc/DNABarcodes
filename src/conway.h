// conway.h

#ifndef __CONWAY_H_INCLUDED__ 
#define __CONWAY_H_INCLUDED__

#include <vector>

#include "distance.h"
#include "sequence.h"

#include <boost/shared_ptr.hpp>

class Conway { 

  private:

  public:
    static std::vector<Sequence> close(const std::vector<Sequence> &, const std::vector<Sequence> &, boost::shared_ptr<Distance>, const unsigned int, const size_t n);
    static std::vector<Sequence> close(const std::vector<Sequence> &, const std::vector<Sequence> &, boost::shared_ptr<Distance>, const unsigned int, const size_t n, const std::vector<Sequence> &);

};

#endif 
