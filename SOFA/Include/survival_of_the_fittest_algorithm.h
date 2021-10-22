#ifndef SOFA_INCLUDE_SURVIAVAL_OF_THE_FITTEST_ALGORITHM_H
#define SOFA_INCLUDE_SURVIAVAL_OF_THE_FITTEST_ALGORITHM_H

#include <cmath>
#include <list>
#include <valarray>
#include <vector>

#include "individual_SOFA.h"
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
  const double dispersion_a;
  const double dispersion_b;
  bool maximize;

  //problem parameters
  unsigned int dimension_;
  std::valarray<Boundaries> parallelepiped_;
  double (*objective_function_)(std::valarray<double>);

  //
  double denominator_;
  double lowest_fitness_;
  IndividualSOFA fittest_;
  std::vector<IndividualSOFA> population_;

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
                     unsigned int dim, unsigned int ipz, double scbd,
                     double da, double db);
  void SetBoundaries(std::valarray<Boundaries>);
  void SetBoundaries(int _size, ...);
  void SetMaximize();
  void SetMinimize();

  std::valarray<double> Optimize();
};
}
#endif //SOFA_INCLUDE_SURVIAVAL_OF_THE_FITTEST_ALGORITHM_H
