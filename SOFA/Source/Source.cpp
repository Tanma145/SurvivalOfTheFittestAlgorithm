#include <cmath>
#include <iomanip>
#include <iostream>
#include <valarray>

#include "../Include/global_optimization.h"
#include "../Include/individual.h"
#include "../Include/test_functions.h"

int main() {
  //set objective function
  SurvivalOfTheFittestAlgotithm::GlobalOptimization sofa(TestFunction2, 2);

  //set parallelepiped
  std::valarray<SurvivalOfTheFittestAlgotithm::Boundaries> bounds = { {0, 10}, {0, 10} };
  sofa.SetBoundaries(bounds);

  sofa.GenerateInitialPopulation(100);

  //output
  std::cout << "initial population" << std::endl;
  for (SurvivalOfTheFittestAlgotithm::Individual cur : sofa.population_) {
      std::cout << std::fixed << std::setprecision(14) << "(" << cur.genotype[0] << ", " << cur.genotype[1] << ") " << cur.fitness << " " << cur.probability << " " << sofa.denominator << std::endl;
  }
  std::cout << std::endl;
  
  for (int i = 0; i < 5000; i++) {
    sofa.GenerateChild();
    std::cout << "(" << sofa.population_.back().genotype[0] << ", " << sofa.population_.back().genotype[1] << ") " << sofa.population_.back().fitness << " " << sofa.population_.back().probability << " " << sofa.denominator << std::endl;
  }

  std::cout << std::endl << "Fittest: " << "(" << sofa.GetFittestGene(0) << ", " << sofa.GetFittestGene(1) << ") " << sofa.GetFittestFitness();
  return 0;
}