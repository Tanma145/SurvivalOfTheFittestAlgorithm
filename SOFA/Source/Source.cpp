#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <valarray>

#include "../Include/survival_of_the_fittest_algorithm.h"
#include "../Include/individual.h"
#include "../Include/test_functions.h"

//#define TEST
#ifdef TEST
int main() {
  std::vector<double> vec = {0, 2, 4, 6, 8, 10, 12, 14, 16 };
  double key;
  while (true){
    std::cout << "Search: ";
    std::cin >> key;
    int left = 0;
    int right = vec.size() - 1;
    int middle;
    while (left != right) {
      middle = floor((left + right) / 2.0);
      if (vec[middle] < key) {
        left = middle + 1;
      }
      else {
        right = middle;
      }
    }
    //if (vec[left] = key) 
    std::cout << "Found: " << vec[right] << std::endl;
  }
  return 0;
}
#else
int main() {
  double initial_population_size = 100;
  double stopping_criteria = 0.005;
  int number_of_tests = 10;
  double (*func)(std::valarray<double>) = TestFunction1;

  std::cout << "number_of_tests = " << number_of_tests << std::endl;
  std::cout << "initial_population_size = " << initial_population_size << std::endl;
  std::cout << "stopping_criteria = " << stopping_criteria << std::endl;
  
  double total_runtime = 0;
  for (int i = 0; i < number_of_tests; i++){
    //calcumations
    int time_start = clock(); //timer start
    GlobalOptimization::SurvivalOfTheFittestAlgorithm sofa(func, 2, initial_population_size, stopping_criteria);
    sofa.SetBoundaries(2, 0.0, 10.0, 0.0, 10.0); //Boundaries must be specified as double (i. e. 2.0 not 2)
    std::valarray<double> fit = sofa.FindMaximum();
    int time_end = clock(); //timer end
    
    //runtime
    double runtime = (time_end - time_start) / (double) CLOCKS_PER_SEC;
    total_runtime += runtime;

    //output
    std::cout.setf(std::ios::fixed);
    std::cout << std::setprecision(8);
    std::cout << "(" << fit[0] << ", " << fit[1] << ")";
    std::cout << " Value: " << func(fit);
    std::cout.unsetf(std::ios::fixed);
    std::cout << "; Population size: " << sofa.population_size_;
    std::cout << "; Runtime: " << runtime << std::endl;
  }
  std::cout << std::endl << "Average runtime (sec): " << total_runtime /  number_of_tests;

  return 0;
}
#endif