#include <cmath>
#include <iomanip>
#include <iostream>
#include <valarray>

#include "../Include/global_optimization.h"
#include "../Include/individual.h"
#include "../Include/test_functions.h"

int main() {
  GlobalOptimization::SurvivalOfTheFittestAlgorithm sofa(TestFunction2, 2, 100, 0.01);
  sofa.SetBoundaries(2, 0.0, 10.0, 0.0, 10.0); //Boundaries must be specified as double (i. e. 2.0 not 2)

  std::cout << std::fixed << std::setprecision(8);
  for (int i = 0; i < 100; i++){
    std::valarray<double> fit = sofa.FindMaximum(); 
    std::cout << "(" << fit[0] << ", " << fit[1] << ")" << std::endl;
  }
  /*
  //set parallelepiped
  std::valarray<GlobalOptimization::Boundaries> bounds = { {0, 10}, {0, 10} };
  sofa.SetBoundaries(bounds);

  sofa.GenerateInitialPopulation(100);

  //output
  std::cout << "initial population" << std::endl;
  for (GlobalOptimization::Individual cur : sofa.population_) {
      std::cout << std::fixed << std::setprecision(14) << "(" << cur.genotype[0] << ", " << cur.genotype[1] << ") " << cur.fitness << " " << cur.probability << " " << sofa.denominator << std::endl;
  }
  std::cout << std::endl;
  
  for (int i = 0; i < 4000; i++) {
    sofa.GenerateChild();
    std::cout << "(" << sofa.population_.back().genotype[0] << ", " << sofa.population_.back().genotype[1] << ") " << sofa.population_.back().fitness << " " << sofa.population_.back().probability << " " << sofa.denominator << std::endl;
  }

  std::cout << std::endl << "Fittest: " << "(" << sofa.GetFittestGene(0) << ", " << sofa.GetFittestGene(1) << ") " << sofa.GetFittestFitness();
  */
  return 0;
}