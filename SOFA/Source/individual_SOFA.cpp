#include "..\Include\individual_SOFA.h"

GlobalOptimization::IndividualSOFA& GlobalOptimization::IndividualSOFA::operator=(GlobalOptimization::IndividualSOFA& ind)  {
  genotype = ind.genotype;
  fitness = ind.fitness;
  probability = ind.probability;
  cumulative_probability = ind.cumulative_probability;
  return *this;
}
