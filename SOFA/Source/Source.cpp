#define _CRT_SECURE_NO_WARNINGS

#include <cmath>
#include <ctime>
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <valarray>

#include "../Include/survival_of_the_fittest_algorithm.h"
#include "../Include/individual.h"
#include "../Include/test_functions.h"

//why yes, i use using, how could you tell?
using std::cout;
using std::ofstream;
using std::endl;

char* settime(struct tm* u)
{
    char s[40];
    char* tmp;
    for (int i = 0; i < 40; i++) s[i] = 0;
    int length = strftime(s, 40, "%d.%m.%Y %H-%M-%S", u);
    tmp = (char*)malloc(sizeof(s));
    strcpy(tmp, s);
    return(tmp);
}

int main() {
  //parameters
  unsigned int initial_population_size = 1;
  double stopping_criteria = 0.00001;
  double dispersion_a = 0.7;
  double dispersion_b = 5e-6;
  int number_of_tests = 10;
  int dimensions = 2;
  std::string title = "TestFunction1";
  double (*func)(std::valarray<double>) = TestFunction1;
  std::valarray<GlobalOptimization::Boundaries> bounds = { {9.0, 10.0}, {3.0, 4.0} };
  // = { {-60.0, 10.0}, {-60.0, 10.0}, {-60.0, 10.0}, {-60.0, 10.0}, {-60.0, 10.0} }

  //output 1
  cout << "objective_function = " << title << endl;
  cout << "number_of_tests = " << number_of_tests << endl;
  cout << "initial_population_size = " << initial_population_size << endl;
  cout << "stopping_criteria = " << stopping_criteria << endl;
  cout << "dispersion(k) =  k ^ -(" << dispersion_a << " + " << dispersion_b << " * k)" << endl;
  cout << "domain of objective function = ";
  cout << "[" << bounds[0].min << "; " << bounds[0].max << "]";
  for (int i = 1; i < dimensions; i++) {
    cout << "x[" << bounds[i].min << "; " << bounds[i].max << "]";
  }
  cout << endl << endl;
  //be sure to change the file path
  std::string filename = "C:\\Users\\tanma\\Documents\\SOFA tests\\SOFA test "; 

  char* file_date;
  const time_t timer = time(NULL);
  file_date = settime(localtime(&timer));
  filename += file_date;
  filename += ".txt";
  ofstream fout(filename, std::ios_base::out);

  fout << "number_of_tests = " << number_of_tests << endl;
  fout << "initial_population_size = " << initial_population_size << endl;
  fout << "stopping_criteria = " << stopping_criteria << endl;
  fout << "dispersion(k) =  k ^ -(" << dispersion_a << " + " << dispersion_b << " * k)" << endl << endl;

  //setting things up
  GlobalOptimization::SurvivalOfTheFittestAlgorithm sofa(func, dimensions, initial_population_size, stopping_criteria, dispersion_a, dispersion_b);
  sofa.SetBoundaries(bounds);
  sofa.SetMaximize();

  double total_runtime = 0.0;
  for (int i = 0; i < number_of_tests; i++){
    //calcumations
    int time_start = clock(); //timer start

    std::valarray<double> fit = sofa.Optimize();

    int time_end = clock(); //timer end
    
    //runtime
    double runtime = (time_end - time_start) / (double) CLOCKS_PER_SEC;
    total_runtime += runtime;

    //output 2
    cout.setf(std::ios::fixed);
    cout << std::setprecision(8);
    cout << "(" << fit[0];
    for (int i = 1; i < dimensions; i++) {
      cout << ", " << fit[i];
    }
    cout << ")";
    cout.unsetf(std::ios::fixed);
    cout << " Value: " << func(fit);
    cout << "; Population size: " << sofa.population_size_;
    cout << "; Runtime: " << runtime << endl;

    fout.setf(std::ios::fixed);
    fout << std::setprecision(8);
    fout << "(" << fit[0];
    for (int i = 1; i < dimensions; i++) {
        fout << ", " << fit[i];
    }
    fout << ")";
    fout.unsetf(std::ios::fixed);
    fout << " Value: " << func(fit);
    fout << "; Population size: " << sofa.population_size_;
    fout << "; Runtime: " << runtime << endl;
  }
  cout << std::endl << "Average runtime (sec): " << total_runtime / number_of_tests;
  fout << std::endl << "Average runtime (sec): " << total_runtime / number_of_tests;

  fout.close();
  return 0;
}