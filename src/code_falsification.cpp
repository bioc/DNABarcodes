//in code_falsification.cpp

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(BH)]]
#include <Rcpp.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "chromosome.h"
#include "chromosome.h" 
#include "conway.h"
#include "distance.h"
#include "genericchromosome.h"
#include "greedyevolution.h"
#include "hammingdistance.h"
#include "levenshteindistance.h"
#include "phaseshiftdist.h"
#include "sequence.h"
#include "sequence.h" 
#include "sequencelevenshteindistance.h"
#include "sequencepool.h"

#include <set>
#include <boost/foreach.hpp>

void calc_mutated_closure(std::set<Sequence> &, const size_t, const std::string&);

std::set<Sequence> convert_to_seq_reads(const std::set<Sequence> &, const size_t);
std::set<Sequence> append_sequence(const Sequence &seq, const size_t n);

void calc_mutated_closure(std::set<Sequence> &mutated_closure, const size_t d, const std::string &metric) {

  // Stop if we cannot mutate anymore
  if (d  == 0) {
    return;
  }

  // Set to store all these cool mutated sequences
  std::set<Sequence> all_new_mutations;

  // Go through every sequence and mutate it
  BOOST_FOREACH(const Sequence &seq, mutated_closure) {
    std::set<Sequence> new_shift_mutations;

    // Delete in the beginning
    new_shift_mutations.insert(seq.remove(0));

    // Insert in the beginning
    BOOST_FOREACH(uint64_t ins_base, Sequence::REAL_BASES) {
      new_shift_mutations.insert(seq.insert(0, ins_base));
    }

    // Calc mutated closures on new shift mutations
    calc_mutated_closure(new_shift_mutations, d-1, metric);

    // Add to all mutations
    all_new_mutations.insert(new_shift_mutations.begin(), new_shift_mutations.end());

    // Substitutions
    std::set<Sequence> new_sub_mutations;

    for (size_t idx = 0; idx < seq.length(); idx++) {
      BOOST_FOREACH(uint64_t sub_base, Sequence::REAL_BASES) {
        new_sub_mutations.insert(seq.substitute(idx, sub_base));
      }
    }

    calc_mutated_closure(new_sub_mutations, d-1, metric);
    all_new_mutations.insert(new_sub_mutations.begin(), new_sub_mutations.end());
  }
  mutated_closure.insert(all_new_mutations.begin(),all_new_mutations.end());
}

std::set<Sequence> append_sequence(const Sequence &seq, const size_t n) {
  std::set<Sequence> appended_set;
  
  if (seq.length() == n) {
    appended_set.insert(seq);
  }

  if (seq.length() >= n) {
    return appended_set;
  }

  BOOST_FOREACH(uint64_t base, Sequence::REAL_BASES) {
    Sequence appended_seq = seq.append(base);
    std::set<Sequence> appended_seq_set = append_sequence(appended_seq, n);
    appended_set.insert(appended_seq_set.begin(), appended_seq_set.end());
  }

  return appended_set;
}

std::set<Sequence> convert_to_seq_reads(const std::set<Sequence> &mutated_closure, const size_t n) {

  std::set<Sequence> mutated_closure_sequences;

  BOOST_FOREACH(const Sequence &seq, mutated_closure) {
    if (seq.length() < n) {
      // Append with every possible combination
      std::set<Sequence> appended_set = append_sequence(seq, n);
      mutated_closure_sequences.insert(appended_set.begin(), appended_set.end());
    } else if (seq.length() > n) {
      // Truncate
       mutated_closure_sequences.insert(seq.truncate(n));
    } else {
      // Just the right length
      mutated_closure_sequences.insert(seq);
    }
  }
  return mutated_closure_sequences;
}

// [[Rcpp::export(".code_falsification")]]
Rcpp::DataFrame code_falsification(const std::vector< std::string > str_barcodes, const unsigned long int d, const std::string metric, const unsigned int cores) {

#ifdef _OPENMP
  omp_set_num_threads(cores);
#endif

  Rcpp::CharacterVector res_barcode;
  Rcpp::CharacterVector res_mutation;

  //#pragma omp parallel for schedule(dynamic) default(none) shared(d, barcodes, n_seqs, dist)

  for(size_t i = 0; i < str_barcodes.size(); i++) {
    std::string seq = str_barcodes.at(i);

    std::set<Sequence> mutated_closure;
    mutated_closure.insert(Sequence(seq));

    // Mutate
    calc_mutated_closure(mutated_closure, d, metric);

    // Produce reads
    std::set<Sequence> mutated_closure_sequences = convert_to_seq_reads(mutated_closure, seq.length());

    // Produce results
    BOOST_FOREACH(const Sequence &read, mutated_closure_sequences) {
      res_barcode.push_back(seq);
      res_mutation.push_back(read.asString());
    }
  }

  return Rcpp::DataFrame::create( Named("barcode")= res_barcode, Named("mutation") = res_mutation, Named("stringsAsFactors") = false);
}

