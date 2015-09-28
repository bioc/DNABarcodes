// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#include <ctime>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/foreach.hpp>
#include <boost/pointer_cast.hpp>

#include "helpers.h"
#include "sequence.h"
#include "sequencepool.h"
#include "conway.h"
#include "distance.h"
#include "cachedistance.h"
#include "phaseshiftdist.h"
#include "greedyevolution.h"
#include "genericchromosome.h"
#include "chromosome.h"
#include "maxclique_pattabiraman_heuristic.h"

#include "create_distance_func.h"

using namespace Rcpp;

std::vector<Sequence> create_dnabarcodes_clique(const std::vector<Sequence> &pool, const unsigned int n, const unsigned long int d, boost::shared_ptr< Distance > dist) {
  size_t N = pool.size();

  MaxClique::Graph g(N);

  Rcout << "2) Calculating distance graph ... ";
  Rcout << std::flush;

  int num_edges = 0;

  bool interrupt = false;    

  size_t loc_d = d;
  // Calculate distance graph
#pragma omp parallel for schedule(dynamic) default(none) shared(pool, g, loc_d, N, num_edges, dist, interrupt)
  for (size_t i = 0; i < N; i++) {
    if (interrupt) continue;

    for (size_t j = i+1; j < N; j++) {
      if (interrupt) continue;

#ifdef _OPENMP
      if (omp_get_thread_num() == 0) // only in master thread!
#endif
        if (check_interrupt()) {
          interrupt = true;
          Rprintf("\nInterrupt detected. Will try to stop gracefully.\n");
          continue;
        }

      if (dist->distance(pool[i],pool[j]) >= loc_d) {
#pragma omp critical
        {
          boost::add_edge(i,j,g);
          ++num_edges;
        }
      }
    }
  }

  if (interrupt) return std::vector<Sequence>();

  Rcout << " done " << std::endl;
  Rcout << std::flush;

  /*
     double density = (double) num_edges / (double) (((N * N) - N) / 2.0);
     Rcout << "Densitiy of the Graph is: " << density << std::endl;
     Rcout << "Hint: " << double(N) / (double) (((N * N) - N) / 2.0) << " would be a sparse graph and 1.0 would be a dense graph." << std::endl;
     Rcout << "Avoid dense graphs. Instead, increase the parameter \"d\" to increase the distance between barcodes." << std::endl;
   */

  // Get clique
  Rcout << "3) Calculating clique ... ";
  Rcout << std::flush;
  std::vector<int> clique = MaxCliquePattabiramanHeuristic::static_max_clique(g, N, 0);
  Rcout << " done " << std::endl;
  Rcout << std::flush;

  /*
     Rcout << "Size of max clique: " << clique->size() << std::endl;
   */

  std::vector<Sequence> clique_sequences;

  BOOST_FOREACH(int element, clique) {
    clique_sequences.push_back(pool[element]);
  }

  return clique_sequences;
}

std::vector<Sequence> create_dnabarcodes_conway(const std::vector<Sequence> &pool, const unsigned int n, const unsigned long int d, boost::shared_ptr< Distance > dist) {
  // FIXME pool not necessary
  // Empty seed on which to close
  std::vector<Sequence> seed;

  Rcout << "2) Conway closing... ";
  Rcout << std::flush;
  std::vector<Sequence> sequences = Conway::close(seed, pool, dist, d, n);
  Rcout << " done " << std::endl;
  Rcout << std::flush;
  return sequences;
}

std::vector<Sequence> create_dnabarcodes_sampling(const std::vector<Sequence> &pool, const unsigned int n, const unsigned long int d, boost::shared_ptr< Distance > dist, const size_t iterations) {
  // Similar to conway, but with random seeds
  size_t max_size = 0;
  std::vector<Sequence> max_set;

  Rcout << "2) Sampling ... ";
  Rcout << std::flush;

  // Our random number generator (FIXME might become a problem with parallelism)

  unsigned int n_seeds  = 3;
  bool interrupt = false;
#pragma omp parallel for schedule(dynamic) default(shared)
  for (size_t i = 0; i < iterations; i++) {
    if (interrupt) continue;

#ifdef _OPENMP
    if (omp_get_thread_num() == 0) // only in master thread!
#endif
      if (check_interrupt()) {
        interrupt = true;
        Rprintf("\nInterrupt detected. Will try to stop gracefully.\n");
        continue;
      }

#ifdef _OPENMP
    boost::random::mt19937_64 gen(omp_get_thread_num() * i);
#else
    boost::random::mt19937_64 gen(i);
#endif
    boost::random::uniform_int_distribution<uint64_t> random_distribution(0, pool.size() - 1);

    std::vector<Sequence> seed;

    unsigned int added    = 0;
    unsigned int failures = 0;

    do {
      unsigned int random_number = random_distribution(gen);
      //std::cerr << "Random number is " << random_number << std::endl;
      Sequence seq_to_add = pool.at(random_number);
      if (dist->is_seq_insertable(seed, seq_to_add, n, d)) {
        seed.push_back(seq_to_add);
        added++;
      } else {
        failures++;
      }
    } while ((added < n_seeds) && (failures < 1000));

    std::vector<Sequence> sequences = Conway::close(seed, pool, dist, d, n);

#pragma omp critical
    {
      if (sequences.size() > max_size) {
        max_size = sequences.size();
        max_set = sequences;

        //Rcout << "New max: " << max_size << std::endl;
      }
    }
  }

  if (interrupt) return std::vector<Sequence>();

  Rcout << " done " << std::endl;
  Rcout << std::flush;

  return max_set;
}

std::vector<Sequence> create_dnabarcodes_ashlock(const std::vector<Sequence> &pool, const unsigned int n, const unsigned long int d, boost::shared_ptr< Distance > dist, const size_t iterations, const size_t population) {
  size_t n_seeds        = 3;

  GreedyEvolutionary ge;

  Rcout << "2) Initiating Chromosomes";
  Rcout << std::flush;

  std::vector< boost::shared_ptr<GenericChromosome> > chromosomes(population);

  BOOST_FOREACH(boost::shared_ptr<GenericChromosome> &chromosome_pointer, chromosomes) {
    chromosome_pointer = boost::shared_ptr<GenericChromosome>(new Chromosome(d, dist, pool, n, n_seeds));
  }
  Rcout << " done " << std::endl;
  Rcout << std::flush;

  Rcout << "3) Running Greedy Evolutionary";
  Rcout << std::flush;
  boost::shared_ptr<Chromosome> best_chromosome = boost::dynamic_pointer_cast<Chromosome>(ge.run(iterations, chromosomes));
  Rcout << " done " << std::endl;
  Rcout << std::flush;

  std::vector<Sequence> sequences;

  if (best_chromosome) {
    sequences = best_chromosome->getCode();
  }

  return sequences;
}

// [[Rcpp::export(".create_dnabarcodes")]]
std::vector< std::string > create_dnabarcodes(
    const unsigned int n,
    const unsigned long int d,
    const std::string metric,
    const std::string generation,
    const bool filter_triplets,
    const bool filter_gc,
    const bool filter_self_complementary,
    const std::vector< std::string > str_pool,
    const unsigned int iterations,
    const unsigned int population,
    const unsigned int cores,
    const bool use_cache,
    const unsigned int cost_sub,
    const unsigned int cost_indel
    ) {

#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif

  // FIXME: pool unnecessary for method "conway"

  /*******************
   * 1) Generate pool *
   ********************/
  Rcout << "1) Creating pool ... ";
  Rcout << std::flush;
  std::vector<Sequence> pool;

  if (str_pool.size() > 0) {
    pool = SequencePool::generate(n, str_pool, filter_triplets, filter_gc, filter_self_complementary);
  } else {
    pool = SequencePool::generate(n, filter_triplets, filter_gc, filter_self_complementary);
  }

  size_t N = pool.size();

  Rcout << " of size " << N << std::endl;
  Rcout << std::flush;

  if (N == 0) {
	Rcpp::stop("Unexpected condition occurred: Empty pool generated.");
  }

  /***********************************
   * Get the right distance algorithm *
   ************************************/

  boost::shared_ptr< Distance > dist, tmp_dist;

  tmp_dist = create_distance_func(metric, cost_sub, cost_indel);

  if (use_cache) {
    Rcout << "1b) Creating cache ... " << std::flush;
    dist = boost::shared_ptr<Distance>(new CacheDistance(tmp_dist)); 
    Rcout << " done " << std::endl << std::flush;
  } else {
    dist = tmp_dist;
  }

  std::vector<Sequence> result_sequences;
  std::vector<std::string> str_sequences;

  if (0 == generation.compare("clique")) {
    result_sequences = create_dnabarcodes_clique(pool, n, d, dist);
  } else if (0 == generation.compare("conway")) {
    result_sequences = create_dnabarcodes_conway(pool, n, d, dist);
  } else if (0 == generation.compare("sampling")) {
    result_sequences = create_dnabarcodes_sampling(pool, n, d, dist, iterations);
  } else if (0 == generation.compare("ashlock")) {
    result_sequences = create_dnabarcodes_ashlock(pool, n, d, dist, iterations, population);
  } else {
    Rcpp::stop("Unrecognized generation method given.");
  }

  BOOST_FOREACH(Sequence seq, result_sequences) {
    str_sequences.push_back(seq.asString());
  }

  return str_sequences;
}

// [[Rcpp::export(".test_distance")]]
unsigned int test_distance(std::string str_seq1, std::string str_seq2) {
  Sequence seq1(str_seq1);
  Sequence seq2(str_seq2);
  clock_t start, end;

  start = std::clock();
  unsigned long int d = PhaseshiftDist::static_distance(seq1, seq2);
  end = std::clock();
  Rcout << "t=" << (double)(end-start)/CLOCKS_PER_SEC << std::endl;
  Rcout << std::flush;
  return d;
}

