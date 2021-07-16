#ifndef SOFA_INCLUDE_GLOBAL_OPTIMIZATION_H
#define SOFA_INCLUDE_GLOBAL_OPTIMIZATION_H

#include <cmath>
#include <list>
#include <valarray>
#include <vector>

#include "Individual.h"

namespace GlobalOptimization{
struct Boundaries {
  double min;
  double max;
};

class SurvivalOfTheFittestAlgorithm {
protected:
  //method parameters
  const double stopping_criteria_by_dispersion_;
  const unsigned int initial_population_size_;

  //problem parameters
  unsigned int dimension_;
  std::valarray<Boundaries> parallelepiped_;
  double (*objective_function_)(std::valarray<double>);

  //
  double denominator_;
  double lowest_fitness_;
  Individual fittest_;
  std::vector<Individual> population_;

  //tools
  double Dispersion(unsigned int) const;
  void CalculateProbabilities();
  void GenerateInitialPopulation(int);
  void GenerateChild();
  void AddIndividual(std::valarray<double>);
  std::valarray<double> GetFittest();
  double GetFittestGene(int);
  double GetFittestFitness();
public:
  double population_size_;
  SurvivalOfTheFittestAlgorithm(double (*objective_function)(std::valarray<double>), 
                     unsigned int dim, unsigned int ipz, double scbd);
  void SetBoundaries(std::valarray<Boundaries>);
  void SetBoundaries(int _size, ...);

  std::valarray<double> FindMaximum();
};
}
#endif //SOFA_INCLUDE_GLOBALOPTIMIZATION_H
