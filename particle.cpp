// particle.cpp

#include "particle.hpp"
#include "config.hpp"

particle::~particle(){}

particle::particle(){

  x = randu<vec>(3);
  x(0) *= LX;
  x(1) *= LY;
  x(2) *= LZ;

  q = quaternion();

}

particle::particle(const particle& particle_to_copy){

  x = particle_to_copy.x;

  q = particle_to_copy.q;

}

void particle::update(const vec& position_update, const vec& lie_algebra_update){

  x += position_update;
  q = lie_exp(lie_algebra_update) * q;

  x(0) -= LX*floor(x(0)/LX);
  x(1) -= LY*floor(x(1)/LY);
  x(2) -= LZ*floor(x(2)/LZ);

}

void particle::update(const vec& position_update, const quaternion& quaternion_update){

  x += position_update;
  q = quaternion_update * q;

  x(0) -= LX*floor(x(0)/LX);
  x(1) -= LY*floor(x(1)/LY);
  x(2) -= LZ*floor(x(2)/LZ);

}

void particle::update(const vec& grand_update){

  x += grand_update.subvec(0,2);
  q = lie_exp(grand_update.subvec(3,5)) * q;

  x(0) -= LX*floor(x(0)/LX);
  x(1) -= LY*floor(x(1)/LY);
  x(2) -= LZ*floor(x(2)/LZ);

}

void particle::accept_configuration_from_master_process(const double *const in){

  x = {in[0], in[1], in[2]};

  q = quaternion(in[3], in[4], in[5], in[6]);

}
