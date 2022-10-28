// collection_of_rigid_bodies.hpp

// =============================================================================
// Include guard
#ifndef MY_COLLECTION_OF_RIGID_BODIES_HEADER_INCLUDED
#define MY_COLLECTION_OF_RIGID_BODIES_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies
class fcm_fluid_solver;

// =============================================================================
// Included dependencies
#include <stdlib.h>
#include <armadillo>
#include "linked_list.hpp"
#include "particle.hpp"
using namespace std;
using namespace arma;

class collection_of_rigid_bodies {

public:

  int blobs_per_body;
  int total_num_blobs;
  mat blob_ref_pos;
  mat all_blob_pos;
  mat ref_KTKinv;
  std::vector<particle> bodies;
  std::vector<particle> midstep_bodies;
  std::vector<particle> temp_bodies;
  vec blob_forces;
  vec temp_blob_forces;
  mat U0;
  double body_radius;
  linked_list list;

  ~collection_of_rigid_bodies();
  collection_of_rigid_bodies();

  void initialise(const int myrank, const int totalnodes);
  void fill_all_blob_pos();
  void fill_all_blob_pos_at_midstep();
  void seed_the_bodies(const int myrank);
  bool advance_to_midstep(const vec& V);
  void calculate_forces(const double t);
  double calculate_potential();
  vec multiply_by_K_at_midstep(const vec& in) const;
  vec multiply_by_KT_at_midstep(const vec& in) const;
  bool end_of_step_update(const vec& V, const double nu);
  void reset_to_start_of_step();
  void write_data_to_file(const double t) const;
  void accept_configuration_from_master_process(const int myrank);
  double calculate_divU(fcm_fluid_solver& fcm);
  bool check_for_overlap();
  bool check_for_channel_wall_overlap();

}; // End of class

#endif // MY_COLLECTION_OF_RIGID_BODIES_HEADER_INCLUDED
