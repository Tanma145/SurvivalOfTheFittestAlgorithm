#ifndef SOFA_INCLUDE_TEST_FUNCTIONS_H
#define SOFA_INCLUDE_TEST_FUNCTIONS_H

#include <valarray>

constexpr double PI = 3.14159265358979323846;

//for SOFA functions must be normalized so 0 < f(x) < 1 otherwise denominator of propability increases indefinitely as population size grows

double TestFunction1(std::valarray<double> argument) {
  //(9.8686984535565490, 3.4023723048523760) - maximum found by DE
  if (argument.size() == 2) {
    double x = argument[0];
    double y = argument[1];
    //0.0022 * (130 + f(x))
     return 0.0022 * (130 + pow(sin(x) + 3 * cos(x) + sin(y) + 3 * cos(y), 2) * (x - 0.5 * y));
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
    return 0.0106 * (52 + (sin(x) + 3 * cos(x) + sin(y) + 3 * cos(y)) * (x - 0.5 * y));
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
        sum1 += A * sin(i * PI * (x - 1 / 2)) * sin(j * PI * (y - 1 / 2)) + B * cos(i * PI * (x - 1 / 2)) * cos(j * PI * (y - 1 / 2));
        sum2 += C * sin(i * PI * (x - 1 / 2)) * sin(j * PI * (y - 1 / 2)) + D * cos(i * PI * (x - 1 / 2)) * cos(j * PI * (y - 1 / 2));
      }
    }
    return sqrt(sum1 * sum1 + sum2 * sum2);
  }
  else {
    throw - 1;
  }
}
#endif //SOFA_INCLUDE_TEST_FUNCTIONS_H