// fcm_fluid_solver.hpp

// =============================================================================
// Include guard
#ifndef MY_FCM_FLUID_SOLVER_HEADER_INCLUDED
#define MY_FCM_FLUID_SOLVER_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies
class collection_of_rigid_bodies;

// =============================================================================
// Included dependencies
#include <armadillo>
#include "fft_solver.hpp"
using namespace arma;

class fcm_fluid_solver {

public:

  // Normal arrays as they interface with the FFT solver
  double *ux;
  double *uy;
  double *uz;
  double *fx;
  double *fy;
  double *fz;
  double *sxx;
  double *syy;
  double *szz;
  double *sxy;
  double *sxz;
  double *syz;

  // Armadillo types for everything else
  vec V;
  vec Vtemp;
  vec Vshear;
  vec ux_rand;
  vec uy_rand;
  vec uz_rand;
  vec ux_shear;
  vec uy_shear;
  vec uz_shear;
  vec qx;
  vec qxsq;
  vec qy;
  vec qysq;
  vec qpad;
  vec qpadsq;
  mat gaussx;
  mat gaussy;
  mat gaussz;
  Mat<int> indx;
  Mat<int> indy;
  Mat<int> indz;
  fft_solver fft;
  mat ref_Minv;
  mat ref_N;

  ~fcm_fluid_solver();
  fcm_fluid_solver();

  void initialise(const int total_num_blobs);
  void set_up_gaussians(const mat& Y);
  void produce_random_velocity_field();
  void produce_shear_velocity_field();
  void calculate_shear_blob_velocities();
  void calculate_random_blob_velocities();
  void set_up_force_distribution(const vec& F);
  void produce_deterministic_velocity_field();
  void calculate_deterministic_blob_velocities();
  void produce_reference_matrices(collection_of_rigid_bodies& all_bodies);

}; // End of fcm_fluid_solver class

#endif // MY_FCM_FLUID_SOLVER_HEADER_INCLUDED
