// parallel_random_numbers.cpp

#include <cmath>
#include "parallel_random_numbers.hpp"

// Look up the "Marsaglia polar method" for details.

double gaussian_random(){
  double x1, x2, z1, w;
  w = 10.0;
  while (w > 1.0 || w == 1){
    x1 = 2.0*sprng() - 1.0;
    x2 = 2.0*sprng() - 1.0;
    w = x1*x1 + x2*x2;
  }
  w = sqrt((-2.0 * log(w))/w);
  z1 = w * x1;
  return z1;
}
