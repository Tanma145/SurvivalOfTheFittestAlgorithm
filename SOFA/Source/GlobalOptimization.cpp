#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <ctime>
#include <iterator>
#include <random>
#include <vector>

#include "../Include/GlobalOptimization.h"

//       1         2         3         4         5         6         7         8         9         0
//34567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890

double InverseProbability(double y, double epsilon, double x_0, double x_min, double x_max) {
  if(0 <= y && y <= 1){
    return epsilon * tan(y * atan((x_max - x_0) / epsilon) + (1 - y) * atan((x_min - x_0) / epsilon)) + x_0;
  }
  else {
    throw 0;
  }
}

SurvivalOfTheFittestAlgotithm::GlobalOptimization::GlobalOptimization(double(*_objective_function)(std::valarray<double>), unsigned int n) {
    objective_function = _objective_function;
    dimension_ = n;
    population_size_ = 0;
    std::list<Individual> pop;
    population_ = pop;
}

void SurvivalOfTheFittestAlgotithm::GlobalOptimization::SetBoundaries(std::valarray<Boundaries> bounds) {
    parallelepiped_ = bounds;
}

double SurvivalOfTheFittestAlgotithm::GlobalOptimization::Dispersion(unsigned int n) const{
  return 1.0 - 2 / M_PI * atan(0.1 * (n - 300));
  //return 2.0;
}

void SurvivalOfTheFittestAlgotithm::GlobalOptimization::CalculateProbabilities(){
  //calculating denomerator for probability
  denominator = 0;
  double sum = 0;
  for (Individual cur : population_) {
    denominator += pow(cur.fitness, population_size_);
  }
  //calculating probabilities of choosing particular individual
  for (auto cur = population_.begin(); cur != population_.end(); cur++) {
    double numerator = pow(cur->fitness, population_size_);
    sum += numerator;
    cur->probability = numerator / denominator;
    cur->cumulative_probability = sum / denominator;
  }
}

void SurvivalOfTheFittestAlgotithm::GlobalOptimization::GenerateInitialPopulation(int n){
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
    ind.fitness = objective_function(ind.genotype);
    if (ind.fitness > fittest.fitness) fittest = ind;
    population.push_back(ind);
  }
  population_size_ = n;
  population_ = population;
  CalculateProbabilities();
}

void SurvivalOfTheFittestAlgotithm::GlobalOptimization::GenerateChild(){
  //choosing random individual with probability of J(x_i)^k / (J(x_1)^k + ... + J(x_k)^k)
  std::mt19937 gen(time(0));
  std::uniform_real_distribution<> urd(0, 1);
  double rnd = urd(gen);

  auto it = population_.begin();
  while (rnd > it->cumulative_probability) {
      it++;
  }
  Individual base = *(it);

  /*
  std::list<Individual>::iterator it, first = population_.begin();
  int count, step;
  count = population_size_;
   
  while (count > 0) {
      it = first;
      step = count / 2;
      std::advance(it, step);
      if (!(rnd < it->cumulative_probability)) {
          first = ++it;
          count -= step + 1;
      }
      else
          count = step;
  }
  Individual base = *(first);

  auto nn = std::upper_bound(population_.begin(), population_.end(), rnd, 
    [](double value, const Individual& ind) {
      return value < ind.cumulative_probability;
    });

  Individual base = *(nn);
  */
  //mutation
  for (int i = 0; i < dimension_; i++) {
    rnd = urd(gen);
    base.genotype[i] = InverseProbability(rnd, Dispersion(population_size_ + 1), base.genotype[i], 
      parallelepiped_[i].min, parallelepiped_[i].max);
  }

  base.fitness = objective_function(base.genotype);
  if (base.fitness > fittest.fitness) fittest = base;
  population_.push_back(base);
  population_size_++;
  CalculateProbabilities();
}

double SurvivalOfTheFittestAlgotithm::GlobalOptimization::GetFittestGene(int i){
  return fittest.genotype[i];
}

double SurvivalOfTheFittestAlgotithm::GlobalOptimization::GetFittestFitness(){
  return fittest.fitness;
}
