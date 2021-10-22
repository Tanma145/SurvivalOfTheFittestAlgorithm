#include <stdarg.h>

#include <algorithm>
#include <cmath>
#include <ctime>
#include <iterator>
#include <random>
#include <vector>

#include "../Include/survival_of_the_fittest_algorithm.h"

double InverseProbability(double y, double epsilon, double x_0, double x_min, double x_max) {
  if (0.0 <= y && y <= 1.0) {
    return epsilon * tan(y * atan((x_max - x_0) / epsilon) + (1.0 - y) * atan((x_min - x_0) / epsilon)) + x_0;
  }
  else {
    throw 0;
  }
}
//       1         2         3         4         5         6         7         8
//345678901234567890123456789012345678901234567890123456789012345678901234567890
GlobalOptimization::SurvivalOfTheFittestAlgorithm::SurvivalOfTheFittestAlgorithm(
  double (*objective_function)(std::valarray<double>), unsigned int dim,
  unsigned int ipz = 100, double scbd = 0.01, double da = 0.4,
  double db = 2.5e-6) : initial_population_size_(ipz),
  stopping_criteria_by_dispersion_(scbd), dispersion_a(da), dispersion_b(db) {

  maximize = true;
  objective_function_ = objective_function;
  dimension_ = dim;
  population_size_ = 0;
  std::vector<IndividualSOFA> pop;
  population_ = pop;
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::SetBoundaries(std::valarray<Boundaries> bounds) {
  parallelepiped_ = bounds;
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::SetMaximize() {
  maximize = true;
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::SetMinimize() {
  maximize = false;
}

double GlobalOptimization::SurvivalOfTheFittestAlgorithm::Dispersion(unsigned int k) const {
  /*
  double kk = k - initial_population_size_;
  return pow(kk, -(dispersion_a + dispersion_b * kk));
  */
  double oof;
  if (k - initial_population_size_ < 1000) oof = 1.0;
  else if (k - initial_population_size_ < 5000) oof = 0.0001;
  else oof = 0.00001;
  return oof;
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::CalculateProbabilities() {

  //calculating denomerator for probability
  denominator_ = 0.0;
  double sum = 0.0;
  for (IndividualSOFA cur : population_) {
    denominator_ += pow((cur.fitness - lowest_fitness_) / (fittest_.fitness - lowest_fitness_), population_size_);
  }
  //calculating probabilities of choosing particular individual
  for (IndividualSOFA cur : population_) {
    double numerator = pow((cur.fitness - lowest_fitness_) / (fittest_.fitness - lowest_fitness_), population_size_);
    sum += numerator;
    cur.probability = numerator / denominator_;
    cur.cumulative_probability = sum / denominator_;
  }
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::GenerateInitialPopulation(int n) {
  std::mt19937 gen(time(0));
  std::vector<IndividualSOFA> population;

  //0th iteration
  IndividualSOFA ind;
  ind.genotype.resize(dimension_);
  for (int j = 0; j < dimension_; j++) {
    std::uniform_real_distribution<> urd(parallelepiped_[j].min, parallelepiped_[j].max);
    ind.genotype[j] = urd(gen);
  }
  ind.fitness = (2 * maximize - 1) * objective_function_(ind.genotype);
  population.push_back(ind);
  //
  fittest_ = population.back();
  lowest_fitness_ = population.back().fitness;

  //remaining iterations
  for (int i = 1; i < n; i++) {
    IndividualSOFA ind;
    ind.genotype.resize(dimension_);
    for (int j = 0; j < dimension_; j++) {
      std::uniform_real_distribution<> urd(parallelepiped_[j].min, parallelepiped_[j].max);
      ind.genotype[j] = urd(gen);
    }
    ind.fitness = (2 * maximize - 1) * objective_function_(ind.genotype);
    if (ind.fitness > fittest_.fitness) fittest_ = ind;
    if (ind.fitness < lowest_fitness_) lowest_fitness_ = ind.fitness;
    population.push_back(ind);
  }
  population_size_ = n;
  population_ = population;
  CalculateProbabilities();
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::GenerateChild() {
  //choosing random individual with probability of J(x_i)^k / (J(x_1)^k + ... + J(x_k)^k)
  std::mt19937 gen(time(0));
  std::uniform_real_distribution<> urd(0.0, 1.0);
  double rnd = urd(gen);

  int left = 0;
  int right = population_size_ - 1;
  int middle;
  while (left != right) {
    middle = floor(left * 0.5 + right * 0.5);
    if (population_[middle].cumulative_probability < rnd) {
      left = middle + 1;
    }
    else {
      right = middle;
    }
  }
  IndividualSOFA base = population_[right];

  //mutation
  for (int i = 0; i < dimension_; i++) {
    rnd = urd(gen);
    base.genotype[i] = InverseProbability(rnd, Dispersion(population_size_ + 1), base.genotype[i],
      parallelepiped_[i].min, parallelepiped_[i].max);
  }

  //maybe that wont work with minimization
  base.fitness = (2.0 * maximize - 1.0) * objective_function_(base.genotype);
  if (base.fitness > fittest_.fitness) fittest_ = base;
  if (base.fitness < lowest_fitness_) lowest_fitness_ = base.fitness;
  population_.push_back(base);
  population_size_++;
  CalculateProbabilities();
}

void GlobalOptimization::SurvivalOfTheFittestAlgorithm::AddIndividual(std::valarray<double> genes) {
  IndividualSOFA ind;
  ind.genotype = genes;
  ind.fitness = (2 * maximize - 1) * objective_function_(genes);
  population_.push_back(ind);
  population_size_++;
  CalculateProbabilities();
}

std::valarray<double> GlobalOptimization::SurvivalOfTheFittestAlgorithm::GetFittest() {
  return fittest_.genotype;
}

double GlobalOptimization::SurvivalOfTheFittestAlgorithm::GetFittestGene(int i) {
  return fittest_.genotype[i];
}

double GlobalOptimization::SurvivalOfTheFittestAlgorithm::GetFittestFitness() {
  return fittest_.fitness;
}

std::valarray<double> GlobalOptimization::SurvivalOfTheFittestAlgorithm::Optimize() {
  GenerateInitialPopulation(initial_population_size_);
  while (Dispersion(population_size_) > stopping_criteria_by_dispersion_) {
    GenerateChild();
  }
  return GetFittest();
}
