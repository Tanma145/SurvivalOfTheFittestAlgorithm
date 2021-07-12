#ifndef SOFA_INCLUDE_GLOBALOPTIMIZATION_H
#define SOFA_INCLUDE_GLOBALOPTIMIZATION_H

#include <cmath>
#include <list>
#include <valarray>

#include "Individual.h"

namespace SurvivalOfTheFittestAlgotithm{
struct Boundaries {
  double min;
  double max;
};

class GlobalOptimization {
protected:
  unsigned int dimension_;
  unsigned int population_size_;
  std::valarray<Boundaries> parallelepiped_;
  double (*objective_function)(std::valarray<double>);
  Individual fittest;

  double Dispersion(unsigned int) const;
  void CalculateProbabilities();
public:
  double denominator;
  std::list<Individual> population_;
  GlobalOptimization(double (*objective_function)(std::valarray<double>), unsigned int);
  void SetBoundaries(std::valarray<Boundaries>);
  void GenerateInitialPopulation(int);
  void GenerateChild();

  double GetFittestGene(int);
  double GetFittestFitness();

  std::valarray<double> FindMaximum();
};
}
#endif //SOFA_INCLUDE_GLOBALOPTIMIZATION_H
