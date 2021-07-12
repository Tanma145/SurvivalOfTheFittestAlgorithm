#define _USE_MATH_DEFINES

#include <algorithm>
#include <cmath>
#include <ctime>
#include <iterator>
#include <random>
#include <vector>

#include "../Include/GlobalOptimization.h"

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

double SurvivalOfTheFittestAlgotithm::GlobalOptimization::Dispersion(unsigned int k) const{
  double a = 0.7;
  double b = 2.5e-6;
  return pow(k, -(a + b * k));
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

void SurvivalOfTheFittestAlgotithm::GlobalOptimization::AddIndividual(std::valarray<double> genes){
  Individual ind;
  ind.genotype = genes;
  ind.fitness = objective_function(genes);
  population_.push_back(ind);
  population_size_++;
  CalculateProbabilities();
}

double SurvivalOfTheFittestAlgotithm::GlobalOptimization::GetFittestGene(int i){
  return fittest.genotype[i];
}

double SurvivalOfTheFittestAlgotithm::GlobalOptimization::GetFittestFitness(){
  return fittest.fitness;
}
