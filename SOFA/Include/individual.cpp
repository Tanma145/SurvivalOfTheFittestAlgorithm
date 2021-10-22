#include "Individual.h"

GlobalOptimization::Individual& GlobalOptimization::Individual::operator=(GlobalOptimization::Individual& ind) {
  genotype = ind.genotype;
  fitness = ind.fitness;
  return *this;
}
