// chromosome.h

#ifndef __CHROMOSOME_H_INCLUDED__ 
#define __CHROMOSOME_H_INCLUDED__

#include "distance.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <boost/shared_ptr.hpp>

#include <vector>
//#include <random>
#include <memory>

#include "genericchromosome.h"

#define CHROMOSOME_MAX_FAILURES 1000

std::ostream& stream_sequences(std::ostream& os, std::vector<Sequence> const &);

class Chromosome : public GenericChromosome { 

  private:
    unsigned int                min_dist_;
    boost::shared_ptr<Distance>   dist_func_;
    const std::vector<Sequence> &pool_;
    unsigned int                n_;
    unsigned int                n_seeds_;

    bool                  is_fitness_known_;
    unsigned int          fitness_;
    std::vector<Sequence> seed_;
    std::vector<Sequence> code_;
    boost::random::mt19937_64 gen_;


  public:
    Chromosome(const unsigned int min_dist, boost::shared_ptr<Distance> dist_func, const std::vector<Sequence> &pool, const unsigned int n, const unsigned int n_seeds);

    Chromosome(const Chromosome& other);

    virtual               ~Chromosome();
    virtual unsigned int  fitness();
    virtual void          initialize();
    virtual unsigned int  getPrecalculatedFitness() const;
    virtual void          mutate();
    virtual bool          operator<(const Chromosome&) const;
    virtual bool          isFitnessKnown() const;
    virtual void          print(std::ostream& os) const;
    
    virtual GenericChromosome *clone() const;

    std::vector<Sequence> getCode() const;

};

#endif 
