// rcross.cpp

#include "rcross.hpp"

arma::mat rcross(const arma::mat& x){

  return {{0,x(2),-x(1)}, {-x(2),0,x(0)}, {x(1),-x(0),0}};

}
