#include <stdarg.h>

#include <algorithm>
#include <cmath>
#include <ctime>
#include <iterator>
#include <random>
#include <vector>

#include "../Include/global_optimization.h"

double InverseProbability(double y, double epsilon, double x_0, double x_min, double x_max) {
  if(0 <= y && y <= 1){
    return epsilon * tan(y * atan((x_max - x_0) / epsilon) + (1 - y) * atan((x_min - x_0) / epsilon)) + x_0;
  }
  else {
    throw 0;
  }
}

GlobalOptimization::SurvivalOfTheFittestAlgorithm::SurvivalOfTheFittestAlgorithm(double (*objective_function)(std::valarray<double>),
                                                                      unsigned int dim, unsigned int ipz = 100, double scbd = 0.01)
                                                                      : initial_population_size_(ipz),
                                                                      stopping_criteria_by_dispersion_(scbd) {
  objective_function_ = objective_function;
  dimension_ = dim;
  population_size_ = 0;
  std::list<Individual> pop;
  population_ = pop;
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::SetBoundaries(std::valarray<Boundaries> bounds) {
  parallelepiped_ = bounds;
}

//Boundaries must be specified as double (e. g. 2.0 not 2)
void GlobalOptimization::SurvivalOfTheFittestAlgorithm::SetBoundaries(int _size, ...){ 
  if (_size % 2 == 0) {
    parallelepiped_.resize(_size);
    va_list argList;
    va_start(argList, 2 *_size);
    for (int i = 0; i < _size; i++) {
      parallelepiped_[i].min = va_arg(argList, double);
      parallelepiped_[i].max = va_arg(argList, double);
    }
  }
  else {
    throw -1;
  }
}

double GlobalOptimization::SurvivalOfTheFittestAlgorithm::Dispersion(unsigned int k) const{
  double a = 0.7;
  double b = 2.5e-6;
  return pow(k, -(a + b * k));
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::CalculateProbabilities(){
  
  //calculating denomerator for probability
  denominator = 0;
  double sum = 0;
  for (Individual cur : population_) {
    denominator += pow(cur.fitness / fittest_.fitness, population_size_);
  }
  //calculating probabilities of choosing particular individual
  for (auto cur = population_.begin(); cur != population_.end(); cur++) {
    double numerator = pow(cur->fitness / fittest_.fitness, population_size_);
    sum += numerator;
    cur->probability = numerator / denominator;
    cur->cumulative_probability = sum / denominator;
  }
  
  /*
  for (Individual cur1 : population_) {
    double inv_prob = 0;
    for (Individual cur2 : population_){
      inv_prob += pow(cur2.fitness / cur1.fitness, population_size_);
    }
    cur1.probability = 1 / inv_prob;
  }
  */
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::GenerateInitialPopulation(int n){
  std::mt19937 gen(time(0));
  double rnd;
  std::list<Individual> population;
  for (int i = 0; i < n; i++) {
    Individual ind;
    ind.genotype.resize(dimension_);
    for (int j = 0; j < dimension_; j++) {
      std::uniform_real_distribution<> urd(parallelepiped_[j].min, parallelepiped_[j].max);
      ind.genotype[j] = urd(gen);
    }
    ind.fitness = objective_function_(ind.genotype);
    if (ind.fitness > fittest_.fitness) fittest_ = ind;
    population.push_back(ind);
  }
  population_size_ = n;
  population_ = population;
  CalculateProbabilities();
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::GenerateChild(){
  //choosing random individual with probability of J(x_i)^k / (J(x_1)^k + ... + J(x_k)^k)
  std::mt19937 gen(time(0));
  std::uniform_real_distribution<> urd(0, 1);
  double rnd = urd(gen);

  auto it = population_.begin();
  while (rnd > it->cumulative_probability) {
      it++;
  }
  Individual base = *(it);

  //mutation
  for (int i = 0; i < dimension_; i++) {
    rnd = urd(gen);
    base.genotype[i] = InverseProbability(rnd, Dispersion(population_size_ + 1), base.genotype[i], 
      parallelepiped_[i].min, parallelepiped_[i].max);
  }

  base.fitness = objective_function_(base.genotype);
  if (base.fitness > fittest_.fitness) fittest_ = base;
  population_.push_back(base);
  population_size_++;
  CalculateProbabilities();
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::AddIndividual(std::valarray<double> genes){
  Individual ind;
  ind.genotype = genes;
  ind.fitness = objective_function_(genes);
  population_.push_back(ind);
  population_size_++;
  CalculateProbabilities();
}

std::valarray<double> GlobalOptimization::SurvivalOfTheFittestAlgorithm::GetFittest(){
  return fittest_.genotype;
}

double GlobalOptimization::SurvivalOfTheFittestAlgorithm::GetFittestGene(int i){
  return fittest_.genotype[i];
}

double GlobalOptimization::SurvivalOfTheFittestAlgorithm::GetFittestFitness(){
  return fittest_.fitness;
}

std::valarray<double> GlobalOptimization::SurvivalOfTheFittestAlgorithm::FindMaximum(){
  GenerateInitialPopulation(initial_population_size_);
  while (Dispersion(population_size_) > stopping_criteria_by_dispersion_) {
    GenerateChild();
  }
  return GetFittest();
}
