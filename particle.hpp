// particle.hpp

// =============================================================================
// Include guard
#ifndef MY_PARTICLE_HEADER_INCLUDED
#define MY_PARTICLE_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include <armadillo>
#include "quaternion.hpp"
using namespace arma;

class particle {

public:

  vec3 x;
  quaternion q;

  ~particle();
  particle();
  particle(const particle& particle_to_copy);

  void update(const vec& position_update, const vec& lie_algebra_update);
  void update(const vec& position_update, const quaternion& quaternion_update);
  void update(const vec& grand_update);
  void accept_configuration_from_master_process(const double *const in);

}; // End of particle class

#endif // MY_PARTICLE_HEADER_INCLUDED
