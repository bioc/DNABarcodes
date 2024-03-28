// in greedyevolutionary.cpp

#include <Rcpp.h>

#include "greedyevolution.h"
#include "genericchromosome.h"

#include <algorithm>
#include <memory>
#include <random>

#include <boost/foreach.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "helpers.h"

static bool comparePtrToChromosome(const boost::shared_ptr<GenericChromosome> &a, const boost::shared_ptr<GenericChromosome> &b) { 
  return (a->getPrecalculatedFitness() < b->getPrecalculatedFitness()); 
}

GreedyEvolutionary::GreedyEvolutionary( ) {

}

boost::shared_ptr<GenericChromosome> GreedyEvolutionary::run(const unsigned int n_evolutions, std::vector<boost::shared_ptr< GenericChromosome > > &chromosomes) {
  size_t n_chromosomes = chromosomes.size();
  std::vector<size_t> histogram;
  // random number generator
  std::random_device rd;
  std::mt19937 g(rd());

  // Initialize chromosomes
  BOOST_FOREACH(boost::shared_ptr<GenericChromosome> chromosome, chromosomes) {
    chromosome->initialize();
  }
 
  // Fitness of chromosomes 
  std::vector<unsigned int> fitness(n_chromosomes);

  size_t max_fitness = 0;
  boost::shared_ptr<GenericChromosome> best_chromosome = chromosomes[1]; // first element is also best chromosome

  bool interrupt = false;

  for (unsigned int evolution(0); evolution < n_evolutions; evolution++) {
    double average_fitness = 0.0;

    if (interrupt) continue;

#pragma omp parallel for schedule(dynamic) default(shared)
    for (unsigned int chrom_idx = 0; chrom_idx < n_chromosomes; chrom_idx++) {

      if (interrupt) continue;

// todo reduce the frequency of checks
#ifdef _OPENMP
      if (omp_get_thread_num() == 0) // only in master thread!
#endif
        if (check_interrupt()) {
          interrupt = true;
          Rprintf("\nInterrupt detected. Will try to stop gracefully.\n");
        }

      // Calculate fitness of every chromosome
      fitness[chrom_idx] = (chromosomes[chrom_idx])->fitness();
#pragma omp critical
      {
        if (fitness[chrom_idx] > max_fitness) {

          max_fitness       = fitness[chrom_idx];
          best_chromosome   = chromosomes[chrom_idx];
        }
      }
#pragma omp critical
      average_fitness += fitness[chrom_idx];
    }

    if (interrupt) continue;

    average_fitness = average_fitness / n_chromosomes;

    // Selection/Mutation step

    std::shuffle(chromosomes.begin(), chromosomes.end(), g);
    
    unsigned int max_group  = n_chromosomes /  GREEDYEVOLUTIONARY_GROUP_SIZE  - 1;

    std::vector< boost::shared_ptr<GenericChromosome> > one_group(GREEDYEVOLUTIONARY_GROUP_SIZE);

    for (unsigned int group(0); group <= max_group; group++) {
      for (unsigned int index(0); index < GREEDYEVOLUTIONARY_GROUP_SIZE; index++) {
        one_group[index] = chromosomes[group * GREEDYEVOLUTIONARY_GROUP_SIZE + index];
      }

      std::sort(one_group.rbegin(), one_group.rend(), comparePtrToChromosome);

      for (unsigned int index(0); index < GREEDYEVOLUTIONARY_GROUP_SIZE / 2; index++) {
        one_group[GREEDYEVOLUTIONARY_GROUP_SIZE - index - 1] = boost::shared_ptr<GenericChromosome>(one_group[index]->clone());
        one_group[GREEDYEVOLUTIONARY_GROUP_SIZE - index - 1]->mutate();
      }

      for (unsigned int index(0); index < GREEDYEVOLUTIONARY_GROUP_SIZE; index++) {
        chromosomes[group * GREEDYEVOLUTIONARY_GROUP_SIZE + index] = one_group[index];
      }
    }

  }

  if (interrupt)
    return(boost::shared_ptr<GenericChromosome>());

  return(best_chromosome);
}

