// greedyevolutionary.h

#ifndef __GREEDYEVOLUTIONARY_H_INCLUDED__ 
#define __GREEDYEVOLUTIONARY_H_INCLUDED__

#define GREEDYEVOLUTIONARY_GROUP_SIZE 4

#include "sequence.h"
#include "genericchromosome.h"

#include <vector>
#include <memory>

class GreedyEvolutionary { 

  public:
    GreedyEvolutionary(); 
    boost::shared_ptr<GenericChromosome> run(const unsigned int, std::vector< boost::shared_ptr< GenericChromosome> > &);

};

#endif 
