// in chromosome.cpp

#include "chromosome.h"
#include "distance.h"
#include "conway.h"

#include <sys/time.h>

#include <vector>
#include <algorithm>
#include <ctime>
//#include <random>
//#include <chrono>

#include <boost/foreach.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

/*
 * The copy constructor
 */

Chromosome::Chromosome( const Chromosome& other ) :  
                          min_dist_(other.min_dist_), 
                          dist_func_(other.dist_func_), 
                          pool_(other.pool_), 
                          n_(other.n_), 
                          n_seeds_(other.n_seeds_),
                          is_fitness_known_(other.is_fitness_known_), 
                          fitness_(other.fitness_), 
                          seed_(other.seed_), 
                          code_(other.code_) {

  static  unsigned int copy_class_iterator = 29;

#pragma omp critical 
  {
  copy_class_iterator++;
  //boost::chrono::duration<> time_now = boost::chrono::high_resolution_clock::now().time_since_epoch();
  struct timeval tv;
  gettimeofday(&tv, NULL);
  uint64_t timeSinceEpoch = (uint64_t)(tv.tv_sec) * 1e6 + (uint64_t)(tv.tv_usec);

  //std::cerr << "Copy constructor called " << std::endl;
#ifdef OMP_H
  uint64_t seed = timeSinceEpoch * copy_class_iterator * (omp_get_thread_num()+1);
#else
  uint64_t seed = timeSinceEpoch * copy_class_iterator;
#endif
  gen_.seed(seed);
  }
}

/*
 * The constructor
 */

Chromosome::Chromosome( const unsigned int min_dist, 
                        boost::shared_ptr<Distance> dist_func, 
                        const std::vector<Sequence> &pool, 
                        const unsigned int n,
                        const unsigned int n_seeds) : 
                          min_dist_(min_dist), 
                          dist_func_(dist_func), 
                          pool_(pool), 
                          n_(n), 
                          n_seeds_(n_seeds),
                          is_fitness_known_(false), 
                          fitness_(0) {
  // Initialize random number generator
  static  unsigned int class_iterator = 1;

#pragma omp critical 
  {
  class_iterator++;

  struct timeval tv;
  gettimeofday(&tv, NULL);
  uint64_t timeSinceEpoch = (uint64_t)(tv.tv_sec) * 1e6 + (uint64_t)(tv.tv_usec);

#ifdef OMP_H
  uint64_t seed = timeSinceEpoch * class_iterator * (omp_get_thread_num()+1);
#else
  uint64_t seed = timeSinceEpoch * class_iterator;
#endif
  gen_.seed(seed);
  }

}

// The destructor
Chromosome::~Chromosome() {

}

void Chromosome::initialize() {
  // reset internal states
  is_fitness_known_ = false;
  fitness_          = 0;
  seed_.clear();
  code_.clear();

  if (n_seeds_ > 0) {
    boost::random::uniform_int_distribution<uint64_t> random_distribution(0, pool_.size() - 1);

    unsigned int added    = 0;
    unsigned int failures = 0;

    do {
      unsigned int random_number = random_distribution(gen_);
      //std::cerr << "Random number is " << random_number << std::endl;
      Sequence seq_to_add = pool_.at(random_number);
      if (dist_func_->is_seq_insertable(seed_, seq_to_add, n_, min_dist_)) {
        //std::cerr << "About to add " << seq_to_add << std::endl;
        seed_.push_back(seq_to_add);
        added++;
      } else {
        failures++;
      }
    } while ((added < n_seeds_) && (failures < CHROMOSOME_MAX_FAILURES));
  }
  code_ = seed_;
}

// Overloaded from "GenericChromosome"
unsigned int Chromosome::fitness() {
  if (is_fitness_known_) {
    return fitness_;
  }

  if (dist_func_->min_set_distance(seed_, n_) < min_dist_) {
    is_fitness_known_ = true;
    fitness_          = 0;
    return 0;
  }

  code_             = Conway::close(seed_, pool_, dist_func_, min_dist_, n_);
  fitness_          = code_.size();
  is_fitness_known_ = true;

  return fitness_;
}

void Chromosome::mutate() {
  boost::random::uniform_int_distribution<> dist(0, pool_.size() - 1);

  for (size_t position(0); position < seed_.size(); position++) {
    if (unif_rand() > 0.5) { 
      bool replaced = false;
      do {
        unsigned int replacement = dist(gen_);
        if (std::find(seed_.begin(), seed_.end(), pool_.at(replacement)) == seed_.end()) {
          seed_[position] = pool_.at(replacement);
          replaced=true;
        }
      } while (!replaced);
      is_fitness_known_ = false;
      fitness_          = 0;
      code_.clear();
    }
  }
}

unsigned int  Chromosome::getPrecalculatedFitness() const {
  return fitness_;
}

bool Chromosome::isFitnessKnown() const {
  return is_fitness_known_;
}

bool Chromosome::operator<(const Chromosome& other) const {
  return (fitness_ < other.getPrecalculatedFitness());
}

std::ostream& stream_sequences(std::ostream& os, std::vector<Sequence> const &seqs) {
  for (size_t seq_idx(0); seq_idx < seqs.size(); seq_idx++) {
    os << "\"" << seqs[seq_idx];
    if (seq_idx != (seqs.size() - 1)) {
      os << "\", ";
    } else {
      os << "\"]";
    }
  }
  return os;
}

void  Chromosome::print(std::ostream& os) const {
  os << "{\n";
 
  os << "\"Fitness\":\t";
  os << fitness_;
  os << ",\n";
   
  os << "\"Seed\":\t[";
  stream_sequences(os, seed_);
  os << ",\n";

  os << "\"Code\":\t[";
  stream_sequences(os, code_);

  os << "\n}";
  return;
}

GenericChromosome* Chromosome::clone() const {
  return new Chromosome(*this);
}

std::vector<Sequence> Chromosome::getCode() const {
  return code_;
}

