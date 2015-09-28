// in genericchromosome.cpp

#include "genericchromosome.h"

GenericChromosome::~GenericChromosome() {

}

void GenericChromosome::initialize() {
  return;
}

unsigned int GenericChromosome::fitness() {
  return 0;
}

void GenericChromosome::mutate() {
  return;
}

unsigned int  GenericChromosome::getPrecalculatedFitness() const {
  return fitness_;
}

bool GenericChromosome::isFitnessKnown() const {
  return is_fitness_known_;
}

bool GenericChromosome::operator<(const GenericChromosome& other) const {
  return (fitness_ < other.getPrecalculatedFitness());
}

void  GenericChromosome::print(std::ostream& os __attribute__ ((unused))) const {
  return;
}

GenericChromosome* GenericChromosome::clone() const {
  return new GenericChromosome(*this);
}

