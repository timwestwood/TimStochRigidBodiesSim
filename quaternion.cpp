// quaternion.cpp

#include <cmath>
#include "quaternion.hpp"
#include "rcross.hpp"
#include "config.hpp"

quaternion::~quaternion(){}

quaternion::quaternion(){

  vec4 q = normalise(randn<vec>(4));

  scalar_part = q(0);
  vector_part = q.subvec(1,3);

}

quaternion::quaternion(const double q0, const double q1, const double q2, const double q3){

  scalar_part = q0;
  vector_part = {q1,q2,q3};

}

quaternion::quaternion(const double q0, const vec& q){

  scalar_part = q0;
  vector_part = q;

}

quaternion::quaternion(const quaternion& quaternion_to_copy){

  scalar_part = quaternion_to_copy.scalar_part;
  vector_part = quaternion_to_copy.vector_part;

}

double quaternion::norm() const {

  return sqrt(scalar_part*scalar_part + dot(vector_part,vector_part));

}

void quaternion::normalise_in_place(){

  double temp = sqrt(scalar_part*scalar_part + dot(vector_part,vector_part));

  scalar_part /= temp;
  vector_part /= temp;

}

quaternion quaternion::conj() const {

  return quaternion(scalar_part, -vector_part);

}

mat33 quaternion::matrix() const {

  mat33 R = eye(3,3);

  double temp = 2.0*vector_part(0)*vector_part(0);
  R(1,1) -= temp;
  R(2,2) -= temp;

  temp = 2.0*vector_part(1)*vector_part(1);
  R(0,0) -= temp;
  R(2,2) -= temp;

  temp = 2.0*vector_part(2)*vector_part(2);
  R(0,0) -= temp;
  R(1,1) -= temp;

  temp = 2.0*vector_part(0)*vector_part(1);
  R(1,0) += temp;
  R(0,1) += temp;

  temp = 2.0*vector_part(0)*vector_part(2);
  R(2,0) += temp;
  R(0,2) += temp;

  temp = 2.0*vector_part(1)*vector_part(2);
  R(1,2) += temp;
  R(2,1) += temp;

  temp = 2.0*scalar_part*vector_part(2);
  R(1,0) += temp;
  R(0,1) -= temp;

  temp = 2.0*scalar_part*vector_part(1);
  R(2,0) -= temp;
  R(0,2) += temp;

  temp = 2.0*scalar_part*vector_part(0);
  R(2,1) += temp;
  R(1,2) -= temp;

  return R;

}

vec quaternion::rotate(const vec& v) const {

  return matrix()*v;

}

void quaternion::sqrt_in_place() {

  if (scalar_part < -0.9999999999){

    scalar_part = 0.0;
    vector_part = {0.0, 0.0, 1.0};

  } else {

    scalar_part = sqrt(0.5*(1.0 + scalar_part));

    vector_part /= 2.0*scalar_part;

  }

}

quaternion lie_exp(const vec& u){

  double theta = norm(u);

  if (theta<1e-10){

    return quaternion(1,0,0,0);

  }else{

    double cs = cos(0.5*theta);
    double sn = sin(0.5*theta);

    return quaternion(cs,sn*u(0)/theta,sn*u(1)/theta,sn*u(2)/theta);

  }

}
