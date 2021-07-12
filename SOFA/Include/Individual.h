#ifndef SOFA_INCLUDE_INDIVIDUAL_H
#define SOFA_INCLUDE_INDIVIDUAL_H

#include <valarray>

namespace SurvivalOfTheFittestAlgotithm{
struct Individual {
  std::valarray<double> genotype;
  double fitness;
  double probability;
  double cumulative_probability;
};
}
#endif //SOFA_INCLUDE_INDIVIDUAL_H
