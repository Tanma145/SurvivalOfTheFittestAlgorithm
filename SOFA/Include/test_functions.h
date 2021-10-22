#ifndef SOFA_INCLUDE_TEST_FUNCTIONS_H
#define SOFA_INCLUDE_TEST_FUNCTIONS_H

#include <valarray>

constexpr double MY_PI = 3.14159265358979323846;
namespace PlanctonFitnessParameters {
  constexpr double c = 140;
  constexpr double c_0 = 60;
  constexpr double c_1 = 40;
  constexpr double alpha = 2;
  constexpr double beta = 2.5e-15;
  constexpr double gamma = 333;
  constexpr double delta = 0.01;
  constexpr double epsilon = 0.13;
  constexpr double sigma_1 = 0.25;
  constexpr double sigma_2 = 0.003;
  constexpr double ksi = 5e-19;
}
double PlanktonStrategy(double time, std::valarray<double> argument) {
  int n = argument.size();
  double height = argument[0];
  for (int i = 1; i < n; i++) {
    height += argument[i] * cos(2 * MY_PI * i * (time));
  }
  return height;
}
double FoodFunction(double height) {
  double food = 0;
  if ((-1) * PlanctonFitnessParameters::c < height && height < 0){
    food = PlanctonFitnessParameters::sigma_1 * (tanh(height + PlanctonFitnessParameters::c_1) + 1);
  }
  return food;
}
double EnergyFunction(double time, std::valarray<double> argument) {
  int n = argument.size();
  double height = 0;
  for (int i = 1; i < n; i++) {
    height -= 2 * MY_PI * i * argument[i] * sin(2 * MY_PI * i * (time));
  }
  return height * height;
 
}
double NaturalDeathFunction(double height) {
  double death = PlanctonFitnessParameters::ksi * cosh(height + PlanctonFitnessParameters::c_0);
  return death;
}
double PredatorHeightFunction(double height) {
  double predator = 0;
  if (-PlanctonFitnessParameters::c < height && height < 0) {
    predator = PlanctonFitnessParameters::sigma_2 * (tanh(height + PlanctonFitnessParameters::c_1) + 1);
  }
  return predator;
}
double PredatorTimeFunction(double time) {
  double predator = 0;
  if (0 < time && time < 1) {
    predator = cos(2 * MY_PI * time) - PlanctonFitnessParameters::epsilon * cos(6 * MY_PI * time);
  }
  return predator;
}
double PlanktonResources(double time, std::valarray<double> argument) {
  double resources = 0;
  resources += PlanctonFitnessParameters::alpha * FoodFunction(PlanktonStrategy(time, argument));
  resources -= PlanctonFitnessParameters::beta  * EnergyFunction(time, argument);
  resources -= PlanctonFitnessParameters::gamma * PredatorHeightFunction(PlanktonStrategy(time, argument)) * PredatorTimeFunction(time);
  resources -= PlanctonFitnessParameters::delta * NaturalDeathFunction(PlanktonStrategy(time, argument));
  return resources;
}
double PlanctonFitness(std::valarray<double> argument){
  int integration_steps = 500;
  const double width = 1 / (double) integration_steps;

  double fit = 0; 
  for (int i = 0; i < integration_steps; i++) {
    const double x1 = i * width;
    const double x2 = (i + 1) * width;
    
    fit += (x2 - x1) / 6.0 * (PlanktonResources(x1, argument) 
                              + 4.0 * PlanktonResources(0.5 * (x1 + x2), argument) 
                              + PlanktonResources(x2, argument));
  }
  return fit;
}

double TestFunction1(std::valarray<double> argument) {
  //(9.8686984535565490, 3.4023723048523760) - maximum found by DE on (0; 10)x(0; 10)
  if (argument.size() == 2) {
    double x = argument[0];
    double y = argument[1];
    //0.0022 * (130 + f(x))
     return (pow(sin(x) + 3 * cos(x) + sin(y) + 3 * cos(y), 2) * (x - 0.5 * y));
  }
  else {
    throw 0;
  }
}
double TestFunction2(std::valarray<double> argument) {
  if (argument.size() == 2) {
    double x = argument[0];
    double y = argument[1];
    //0.0106 * (52 + f(x))
    return (52 + (sin(x) + 3 * cos(x) + sin(y) + 3 * cos(y)) * (x - 0.5 * y));
  }
  else {
    throw 0;
  }
}
double RosenbrockFunction(std::valarray<double> argument) {
  if (argument.size() == 2) {
    double x = argument[0];
    double y = argument[1];
    return pow((1 - x), 2) + 100 * pow((y - x * x), 2);
  }
  else{
    throw - 1;
  }
}
double IcicleFunction(std::valarray<double> argument) {
  if (argument.size() == 2) {
    double x = argument[0];
    double y = argument[1];
    return (1 + sin(10 * x) + cos(2 * x) + cos(2 * x + 2 * y) + cos(2 * y) + sin(20 * y) + y * y);
  }
  else{
    throw - 1;
  }
}
double HorrificFunction(std::valarray<double> argument) {
  if (argument.size() == 2) {
    double x = argument[0];
    double y = argument[1];
    double result = 0;
    double sum1 = 0, sum2 = 0;
    double A, B, C, D;
    std::valarray<std::valarray<double>> A_matrix = {{-0.940, -0.536, -0.743},
													 {-0.502,  0.804,  0.769},
													 {-0.428, -0.789,  0.204}};
    std::valarray<std::valarray<double>> B_matrix = {{ 0.590,  0.160, -0.681},
													 { 0.387,  0.945, -0.195},
													 {-0.231,  0.152,  0.295}};
	std::valarray<std::valarray<double>> C_matrix = {{-0.896, -0.613, -0.463},
													 { 0.038, -0.428, -0.714},
													 { 0.103,  0.741, -0.317}};
    std::valarray<std::valarray<double>> D_matrix = {{-0.754, -0.558, -0.989},
													 {-0.702,  0.881,  0.397},
													 {-0.056,  0.085, -0.616}};
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        A = A_matrix[i][j];
        B = B_matrix[i][j];
        C = C_matrix[i][j];
        D = D_matrix[i][j];
        sum1 += A * sin(i * MY_PI * (x - 1 / 2)) * sin(j * MY_PI * (y - 1 / 2)) + B * cos(i * MY_PI * (x - 1 / 2)) * cos(j * MY_PI * (y - 1 / 2));
        sum2 += C * sin(i * MY_PI * (x - 1 / 2)) * sin(j * MY_PI * (y - 1 / 2)) + D * cos(i * MY_PI * (x - 1 / 2)) * cos(j * MY_PI * (y - 1 / 2));
      }
    }
    return sqrt(sum1 * sum1 + sum2 * sum2);
  }
  else {
    throw - 1;
  }
}
#endif //SOFA_INCLUDE_TEST_FUNCTIONS_H