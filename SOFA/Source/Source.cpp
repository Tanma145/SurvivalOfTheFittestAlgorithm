#include <cmath>
#include <iomanip>
#include <iostream>
#include <valarray>

#include "../Include/GlobalOptimization.h"
#include "../Include/Individual.h"

double TestFunction1(std::valarray<double> parameters) {
  //(9.8686984535565490, 3.4023723048523760) - maximum found by DE
  if (parameters.size() == 2) {
    double x = parameters[0];
    double y = parameters[1];
    //function must be normalized so 0 < f(x) < 1 otherwise denominator increases indefinitely as population size grows
    //0.0022 * (130 + f(x))
    return 0.0022 * (130 + pow(sin(x) + 3 * cos(x) + sin(y) + 3 * cos(y), 2) * (x - 0.5 * y));
  }
  else {
    throw 0;
  }
}
double TestFunction2(std::valarray<double> parameters) {
  if (parameters.size() == 2) {
    double x = parameters[0];
    double y = parameters[1];
    return 0.009 * (60 + (sin(x) + 3 * cos(x) + sin(y) + 3 * cos(y)) * (x - 0.5 * y));
  }
  else {
    throw 0;
  }
}

int main() {
  //set objective function
  SurvivalOfTheFittestAlgotithm::GlobalOptimization sofa(TestFunction2, 2);

  //set parallelepiped
  std::valarray<SurvivalOfTheFittestAlgotithm::Boundaries> bounds = { {0, 10}, {0, 10} };
  sofa.SetBoundaries(bounds);

  sofa.GenerateInitialPopulation(10);

  //output
  std::cout << "initial population" << std::endl;
  for (SurvivalOfTheFittestAlgotithm::Individual cur : sofa.population_) {
      std::cout << std::fixed << std::setprecision(14) << "(" << cur.genotype[0] << ", " << cur.genotype[1] << ") " << cur.fitness << " " << cur.probability << " " << sofa.denominator << std::endl;
  }
  std::cout << std::endl;
  
  for (int i = 0; i < 2000; i++) {
    sofa.GenerateChild();
    std::cout << "(" << sofa.population_.back().genotype[0] << ", " << sofa.population_.back().genotype[1] << ") " << sofa.population_.back().fitness << " " << sofa.population_.back().probability << " " << sofa.denominator << std::endl;
  }

  std::cout << std::endl << "Fittest: " << "(" << sofa.GetFittestGene(0) << ", " << sofa.GetFittestGene(1) << ") " << sofa.GetFittestFitness();
  return 0;
}