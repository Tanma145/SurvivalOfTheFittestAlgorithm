#ifndef SOFA_INCLUDE_INDIVIDUAL_H
#define SOFA_INCLUDE_INDIVIDUAL_H

#include <valarray>

namespace GlobalOptimization{
class Individual {
public:
  std::valarray<double> genotype;
  double fitness;
  Individual& operator=(Individual& ind);
};
}
#endif //SOFA_INCLUDE_INDIVIDUAL_H
