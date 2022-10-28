// gmres_solver.hpp

// =============================================================================
// Include guard
#ifndef MY_GMRES_SOLVER_HEADER_INCLUDED
#define MY_GMRES_SOLVER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies
class fcm_fluid_solver;
class collection_of_rigid_bodies;

// =============================================================================
// Included dependencies
#include <armadillo>
using namespace arma;

class gmres_solver {

public:

  int num_blobs;
  int system_size;
  vec RHS;
  mat Q;
  vec beta;
  mat H;
  vec SN;
  vec CS;

  ~gmres_solver();
  gmres_solver();

  void initialise(const int total_num_blobs);
  void assemble_RHS(const collection_of_rigid_bodies& all_bodies, const fcm_fluid_solver& fcm);
  vec multiply_by_Pinv(const vec& in, const collection_of_rigid_bodies& all_bodies, const fcm_fluid_solver& fcm);
  vec multiply_by_saddle_point_matrix(const vec& in, const collection_of_rigid_bodies& all_bodies, fcm_fluid_solver& fcm);
  vec invert_saddle_point_system(const collection_of_rigid_bodies& all_bodies, fcm_fluid_solver& fcm, const int myrank, const double t);

}; // End of class

#endif // MY_GMRES_SOLVER_HEADER_INCLUDED
