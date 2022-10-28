// collection_of_rigid_bodies.cpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <mpi.h>
#include "fcm_fluid_solver.hpp"
#include "collection_of_rigid_bodies.hpp"
#include "rcross.hpp"
#include "config.hpp"

collection_of_rigid_bodies::~collection_of_rigid_bodies(){}
collection_of_rigid_bodies::collection_of_rigid_bodies(){}

void collection_of_rigid_bodies::initialise(const int myrank, const int totalnodes){

  #if SIMULATING_BLOBS

    blobs_per_body = 1;

  #elif SIMULATING_RODS

    blobs_per_body = 22; // Should be twice an odd number for axi-symmetry.

  #elif SIMULATING_TETRAHEDRA

    blobs_per_body = 4;

  #elif SIMULATING_CZECH_HEDGEHOGS

    blobs_per_body = 19;

  #endif

  blob_ref_pos = mat(3,blobs_per_body,fill::zeros);

  total_num_blobs = blobs_per_body*NUM_BODIES;

  list.initialise(total_num_blobs, myrank, totalnodes);

  all_blob_pos = mat(total_num_blobs, 3);

  #if (!MARKOV_CHAIN_MONTE_CARLO || IC_SEED_RANDOM)

  blob_forces = vec(3*total_num_blobs);
  temp_blob_forces = vec(3*total_num_blobs);

  #endif

  #if !MARKOV_CHAIN_MONTE_CARLO

  U0 = mat(6,NUM_BODIES);

  #endif

  #if SIMULATING_RODS

    int pair = 0;

    for (int blob_id=0; blob_id<blobs_per_body; blob_id+=2){

      const double theta = pair*0.5*PI;
      const double x = sqrt(2.0)*BLOB_RADIUS*(pair - 0.25*blobs_per_body + 0.5);

      blob_ref_pos.col(blob_id) = vec({x, BLOB_RADIUS*cos(theta), BLOB_RADIUS*sin(theta)});
      blob_ref_pos.col(blob_id + 1) = vec({x, -BLOB_RADIUS*cos(theta), -BLOB_RADIUS*sin(theta)});

      pair++;

    }

  #elif SIMULATING_TETRAHEDRA

    double fac = 1.05*BLOB_RADIUS;

    blob_ref_pos.col(0) = vec({fac,0,-fac/sqrt(2.0)});
    blob_ref_pos.col(1) = vec({-fac,0,-fac/sqrt(2.0)});
    blob_ref_pos.col(2) = vec({0,fac,fac/sqrt(2.0)});
    blob_ref_pos.col(3) = vec({0,-fac,fac/sqrt(2.0)});

  #elif SIMULATING_CZECH_HEDGEHOGS

    const double d = 2.0*BLOB_RADIUS;

    blob_ref_pos(0, 1) = d;
    blob_ref_pos(0, 2) = 2.0*d;
    blob_ref_pos(0, 3) = 3.0*d;
    blob_ref_pos(0, 4) = -d;
    blob_ref_pos(0, 5) = -2.0*d;
    blob_ref_pos(0, 6) = -3.0*d;

    double theta = 0.0;

    for (int blob_id=7; blob_id < blobs_per_body; blob_id++){

      const int fac = 1 + ((blob_id - 7)/4);

      blob_ref_pos(1, blob_id) = fac*d*cos(theta);
      blob_ref_pos(2, blob_id) = fac*d*sin(theta);

      theta += 0.5*PI;

    }

  #endif

  bodies = std::vector<particle>(NUM_BODIES);
  temp_bodies = bodies;

  body_radius = 0.0;
  for (int i=0; i<blobs_per_body; i++){

    body_radius = std::max<double>(body_radius, norm(blob_ref_pos.col(i)));

  }

  #if !SIMULATING_BLOBS

    mat K(3*blobs_per_body, 6);

    for (int m=0; m < blobs_per_body; m++){

      K(3*m, 0, size(3,3)) = eye(3,3);
      K(3*m, 3, size(3,3)) = rcross(blob_ref_pos.col(m));

    }

    ref_KTKinv = inv(K.t() * K);

  #endif

  #if !MARKOV_CHAIN_MONTE_CARLO

    midstep_bodies = bodies;

    /* Change the position and orientation of the first body for the single-body mobility calculation.
    This change will be overwritten in the call to seed_the_bodies. */

    #if SLIP_CHANNEL

      bodies[0].x = {0.5*LX, 0.25*LY, 0.5*LZ};

    #else

      bodies[0].x = {0.5*LX, 0.5*LY, 0.5*LZ};

    #endif

    bodies[0].q = quaternion(1.0, 0.0, 0.0, 0.0);

  #endif

}

void collection_of_rigid_bodies::fill_all_blob_pos(){

  for (int n=0; n < NUM_BODIES; n++){

    const mat33 R = bodies[n].q.matrix();

    for (int m=0; m < blobs_per_body; m++){

      const vec3 p = bodies[n].x + R*blob_ref_pos.col(m);

      all_blob_pos(n*blobs_per_body + m, 0) = p(0) - LX*floor(p(0)/LX);
      all_blob_pos(n*blobs_per_body + m, 1) = p(1) - LY*floor(p(1)/LY);
      all_blob_pos(n*blobs_per_body + m, 2) = p(2) - LZ*floor(p(2)/LZ);

    }

  }

}

void collection_of_rigid_bodies::fill_all_blob_pos_at_midstep(){

  for (int n=0; n < NUM_BODIES; n++){

    const mat33 R = midstep_bodies[n].q.matrix();

    for (int m=0; m < blobs_per_body; m++){

      const vec3 p = midstep_bodies[n].x + R*blob_ref_pos.col(m);

      all_blob_pos(n*blobs_per_body + m, 0) = p(0) - LX*floor(p(0)/LX);
      all_blob_pos(n*blobs_per_body + m, 1) = p(1) - LY*floor(p(1)/LY);
      all_blob_pos(n*blobs_per_body + m, 2) = p(2) - LZ*floor(p(2)/LZ);

    }

  }

}

bool collection_of_rigid_bodies::check_for_overlap(){

  bool temp = list.check_for_overlap(all_blob_pos); // Returns true if it finds any overlap.

  bool out;

  MPI_Allreduce(&temp, &out, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

  return out;

}

bool collection_of_rigid_bodies::check_for_channel_wall_overlap(){
// TO-DO: Change this so that either the work is shared across processes or the master process just does it and broadcasts.
  #if SLIP_CHANNEL

  bool temp = false;

  for (int n = 0; n < total_num_blobs; n++){

    const double y = all_blob_pos(n, 1);

    if ((y == 0.0) || (y >= 0.5*LY)){

      temp = true;

      break;

    }

  }

  bool out;

  MPI_Allreduce(&temp, &out, 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

  return out;

  #else

  return false;

  #endif

}

void collection_of_rigid_bodies::seed_the_bodies(const int myrank){

  if (myrank == 0){

    ofstream output_data_file(SIMULATION_NAME+std::string("_blob_references.dat"));
    output_data_file << std::scientific << std::setprecision(6);

    for (int i=0; i<blobs_per_body; i++){
      output_data_file << blob_ref_pos(0,i) << " ";
    }
    output_data_file << endl;
    for (int i=0; i<blobs_per_body; i++){
      output_data_file << blob_ref_pos(1,i) << " ";
    }
    output_data_file << endl;
    for (int i=0; i<blobs_per_body; i++){
      output_data_file << blob_ref_pos(2,i) << " ";
    }
    output_data_file << endl;

    output_data_file.close();

  }

  #if READ_IN_INITIAL_CONDITIONS

  if (myrank == 0){

    cout << "Reading in initial conditions..." << endl << endl;

    ifstream state_file(std::string(INITIAL_CONDITIONS_FILE_NAME));

    // Work out how many states are in the file
    double temp;
    int n = 0;

    while (state_file >> temp){

      for (int i=0; i<NUM_BODIES; i++){

        double x, y, z, q0, q1, q2, q3;

        state_file >> x >> y >> z >> q0 >> q1 >> q2 >> q3;

      }

      n++;

    }

    state_file.close();
    state_file.clear();

    #if IC_READ_IN_RANDOM

    const int sample_to_use = randi<int>(distr_param(1, n));

    cout << "Input file contains " << n << " state(s). State " << sample_to_use << " has been selected." << endl << endl;

    #elif IC_READ_IN_LAST

    const int sample_to_use = n;

    #endif

    state_file.open(std::string(INITIAL_CONDITIONS_FILE_NAME));

    for (int m = 0; m < sample_to_use; m++){

      state_file >> temp;

      for (int i=0; i<NUM_BODIES; i++){

        double x, y, z, q0, q1, q2, q3;

        state_file >> x >> y >> z >> q0 >> q1 >> q2 >> q3;

        if (m == sample_to_use-1){

          bodies[i].x = {x, y, z};
          bodies[i].q = quaternion(q0, q1, q2, q3);

        }

      }

      if (m == sample_to_use-1){

        break;

      }

    }

    state_file.close();

  }

  accept_configuration_from_master_process(myrank);

  #else

  if (myrank == 0){

    cout << "Forming initial conditions..." << endl << endl;

    // Return the first body to a random state if we moved to form reference matrices.
    #if !MARKOV_CHAIN_MONTE_CARLO

      bodies[0] = particle();

    #endif

    #if SLIP_CHANNEL

      #if SIMULATING_BLOBS

        quaternion q(1.0, 0.0, 0.0, 0.0);

      #else

        const int domain_dir_max = ((LX >= 0.5*LY) && (LX >= LZ)) ? 0 : (((0.5*LY >= LX) && (0.5*LY >= LZ)) ? 1 : 2);
        const int domain_dir_min = ((LX <= 0.5*LY) && (LX <= LZ)) ? 0 : (((0.5*LY <= LX) && (0.5*LY <= LZ)) ? 1 : 2);

        double lx = max(blob_ref_pos.row(0)) - min(blob_ref_pos.row(0));
        double ly = max(blob_ref_pos.row(1)) - min(blob_ref_pos.row(1));
        double lz = max(blob_ref_pos.row(2)) - min(blob_ref_pos.row(2));

        const int body_dir_max = ((lx >= ly) && (lx >= lz)) ? 0 : (((ly >= lx) && (ly >= lz)) ? 1 : 2);
        const int body_dir_min = ((lx <= ly) && (lx <= lz)) ? 0 : (((ly <= lx) && (ly <= lz)) ? 1 : 2);

        vec3 v1 = {0.0, 0.0, 0.0};
        v1(body_dir_max) = 1.0;
        vec3 v2 = {0.0, 0.0, 0.0};
        v2(domain_dir_max) = 1.0;

        quaternion q(dot(v1, v2), cross(v1, v2));
        q.sqrt_in_place();

        v1 = (q.matrix()).col(body_dir_min);

        v2.zeros();
        v2(domain_dir_min) = 1.0;

        quaternion qtemp(dot(v1, v2), cross(v1, v2));
        qtemp.sqrt_in_place();
        q = qtemp * q;

      #endif

      const mat temp = q.matrix()*blob_ref_pos;
      const double y_lower_bound = std::min(0.25*LY, 0.0 - 1.1*min(temp.row(1))); // Minimum ought to be negative.
      const double y_upper_bound = std::max(0.25*LY, 0.5*LY - 1.1*max(temp.row(1))); // Maximum ought to be positive.

      for (int n = 0; n < NUM_BODIES; n++){

        bodies[n].x(1) -= LY*floor(2.0*bodies[n].x(1)/LY);

      }

      fill_all_blob_pos();

      for (int n = 0; n < NUM_BODIES; n++){

        for (int m = 0; m < blobs_per_body; m++){

          const double y = all_blob_pos(n*blobs_per_body + m, 1);

          if ((y == 0.0) || (y >= 0.5*LY)){

            bodies[n].q = q;
            bodies[n].x(1) = y_lower_bound + (y_upper_bound - y_lower_bound)*randu<double>();

            break;

          }

        }

      }

    #endif

  }

  accept_configuration_from_master_process(myrank);

  fill_all_blob_pos();
  list.populate(all_blob_pos);
  int num_attempts = 0;

  // Assuming the above worked, no blob centres are beyond the walls (if they exist).
  bool bad_setup = check_for_overlap();

  const double force_fac = 1.0/(6.0*PI*FLUID_VISCOSITY*BLOB_RADIUS);

  const double dt = 0.1*BLOB_DIFFUSIVE_TIME;

  while (bad_setup && (num_attempts < 1000)){

    temp_bodies = bodies;

    temp_blob_forces.zeros();

    list.apply_barrier_forces(temp_blob_forces, all_blob_pos);

    MPI_Allreduce(temp_blob_forces.memptr(), blob_forces.memptr(), 3*total_num_blobs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    #if SLIP_CHANNEL

    // This force should stop anything getting too close to the channel walls
    const double L = 2.2*BLOB_RADIUS;
    const double K = (40.0*ENERGY_SCALE/BLOB_RADIUS)/L;

    for (int i=0; i<total_num_blobs; i++){

      const double y = all_blob_pos(i,1);

      if (y <= L){

        blob_forces(3*i + 1) += K*(L-y);

      } else if ((0.5*LY - y <= L) && (y <= 0.5*LY)){

        blob_forces(3*i + 1) -= K*(L - 0.5*LY + y);

      }

    }

    #endif

    #if SIMULATING_BLOBS

      for (int n = 0; n < NUM_BODIES; n++){

        bodies[n].update(force_fac*dt*blob_forces.subvec(3*n, 3*n + 2), vec({0.0, 0.0, 0.0}));

      }

    #else

      for (int n = 0; n < NUM_BODIES; n++){

        vec3 body_force = {0.0, 0.0, 0.0};
        vec3 body_torque = {0.0, 0.0, 0.0};

        const mat33 R = bodies[n].q.matrix();

        for (int m = 0; m < blobs_per_body; m++){

          body_force += R.t() * blob_forces.subvec(3*(n*blobs_per_body + m), 3*(n*blobs_per_body + m) + 2);
          body_torque += cross(blob_ref_pos.col(m), R.t() * blob_forces.subvec(3*(n*blobs_per_body + m), 3*(n*blobs_per_body + m) + 2));

        }

        vec update = force_fac*dt*ref_KTKinv*join_vert(body_force, body_torque);
        update.subvec(0,2) = R*update.subvec(0,2);
        update.subvec(3,5) = R*update.subvec(3,5);

        bodies[n].update(update);

      }

    #endif

    accept_configuration_from_master_process(myrank);

    fill_all_blob_pos();

    #if SLIP_CHANNEL

      for (int n = 0; n < NUM_BODIES; n++){

        bool reset = false;

        for (int m = 0; m < blobs_per_body; m++){

          const double y = all_blob_pos(n*blobs_per_body + m, 1);

          if ((y == 0.0) || (y >= 0.5*LY)){

            reset = true;

            break;

          }

        }

        if (reset){

          bodies[n] = temp_bodies[n];

        }

      }

    #endif

    fill_all_blob_pos();
    list.populate(all_blob_pos);

    bad_setup = check_for_overlap();

    num_attempts++;

    if (myrank == 0){

      cout << "Completed drag-initialisation step " << num_attempts << endl;

    }

  }

  accept_configuration_from_master_process(myrank);

  #endif

  if (myrank == 0){
    cout << "Initial configuration has been distributed by the master process..." << endl << endl;
  }

}

bool collection_of_rigid_bodies::advance_to_midstep(const vec& V){

  midstep_bodies = bodies;

  for (int n=0; n<NUM_BODIES; n++){

    #if SIMULATING_BLOBS

      U0.col(n) = vec({V(3*n), V(3*n + 1), V(3*n + 2), 0.0, 0.0, 0.0});

      midstep_bodies[n].update(0.5*DT*U0.col(n));

    #else

      vec b = vec(6,fill::zeros);

      const mat Q = bodies[n].q.matrix();

      for (int i=0; i<blobs_per_body; i++){

        int id = n*blobs_per_body + i;

        vec3 vel = V.subvec(3*id, 3*id + 2);

        b += join_vert(vel, cross(Q*blob_ref_pos.col(i), vel));

      }

      b.subvec(0,2) = Q.t() * b.subvec(0,2);
      b.subvec(3,5) = Q.t() * b.subvec(3,5);

      U0.col(n) = ref_KTKinv*b;

      U0(0, n, size(3,1)) = Q * U0(0, n, size(3,1));
      U0(3, n, size(3,1)) = Q * U0(3, n, size(3,1));

      midstep_bodies[n].update(0.5*DT*U0.col(n));

    #endif

  }

  fill_all_blob_pos_at_midstep();
  list.populate(all_blob_pos);

  #if SKIP_ON_OVERLAP

  return check_for_channel_wall_overlap() || check_for_overlap();

  #else

  return check_for_channel_wall_overlap();

  #endif

}

void collection_of_rigid_bodies::calculate_forces(const double t){

  #if (EULER_MARUYAMA && !SIMULATING_BLOBS)

    midstep_bodies = bodies; // So the correct values are used in the GMRES solve.

  #endif

  temp_blob_forces.zeros();

  // Start with barrier forces
  list.apply_barrier_forces(temp_blob_forces, all_blob_pos);

  MPI_Allreduce(temp_blob_forces.memptr(), blob_forces.memptr(), 3*total_num_blobs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Apply any other forces

  #if SLIP_CHANNEL

  // This force should stop anything getting too close to the channel walls
  const double L = 2.2*BLOB_RADIUS;
  const double K = (40.0*ENERGY_SCALE/BLOB_RADIUS)/L;

  for (int i=0; i<total_num_blobs; i++){

    const double y = all_blob_pos(i,1);

    if (y <= L){

      blob_forces(3*i + 1) += K*(L-y);

    } else if ((0.5*LY - y <= L) && (y <= 0.5*LY)){

      blob_forces(3*i + 1) -= K*(L - 0.5*LY + y);

    }

  }

  #endif

}

double collection_of_rigid_bodies::calculate_potential(){

  double temp_potential = 0.0;

  // Start with interaction potentials
  list.apply_barrier_potential(temp_potential, all_blob_pos);

  double potential;

  MPI_Allreduce(&temp_potential, &potential, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  // Apply any other potentials

  #if SLIP_CHANNEL

  // This potential should stop anything getting too close to the channel walls
  const double L = 2.2*BLOB_RADIUS;
  const double K = 40.0*ENERGY_SCALE/(L*BLOB_RADIUS);

  for (int i=0; i<total_num_blobs; i++){

    const double y = all_blob_pos(i,1);

    if (y <= L){

      potential += 0.5*K*(L-y)*(L-y);

    } else if ((0.5*LY - y <= L) && (y <= 0.5*LY)){

      potential += 0.5*K*(L - 0.5*LY + y)*(L - 0.5*LY + y);

    }

  }

  #endif

  return potential;

}

vec collection_of_rigid_bodies::multiply_by_K_at_midstep(const vec& in) const {

  const int N = 3*total_num_blobs;

  vec out = vec(N,fill::zeros);

  for (int i=0; i<NUM_BODIES; i++){

    int ipos = 6*i;

    const mat R = midstep_bodies[i].q.matrix();

    for (int j=0; j<blobs_per_body; j++){

      int jpos = 3*(blobs_per_body*i + j);

      out.subvec(jpos,jpos+2) += in.subvec(ipos,ipos+2) + cross(in.subvec(ipos+3,ipos+5), R*blob_ref_pos.col(j));

    }

  }

  return out;

}

vec collection_of_rigid_bodies::multiply_by_KT_at_midstep(const vec& in) const {

  const int N = 6*NUM_BODIES;

  vec out = vec(N,fill::zeros);

  for (int i=0; i<NUM_BODIES; i++){

    int ipos = 6*i;

    const mat R = midstep_bodies[i].q.matrix();

    for (int j=0; j<blobs_per_body; j++){

      int jpos = 3*(blobs_per_body*i + j);

      out.subvec(ipos,ipos+5) += join_vert(in.subvec(jpos,jpos+2), cross(R*blob_ref_pos.col(j), in.subvec(jpos,jpos+2)));

    }

  }

  return out;

}

bool collection_of_rigid_bodies::end_of_step_update(const vec& V, const double nu){

  temp_bodies = bodies;

  for (int n=0; n<NUM_BODIES; n++){

    #if SIMULATING_BLOBS

      bodies[n].update(nu*DT*vec({V(3*n), V(3*n + 1), V(3*n + 2), 0.0, 0.0, 0.0}));

    #else

      bodies[n].update(nu*DT*V.subvec(6*n,6*n+5));

    #endif

  }

  fill_all_blob_pos();
  list.populate(all_blob_pos);

  #if SKIP_ON_OVERLAP

  return check_for_channel_wall_overlap() || check_for_overlap();

  #else

  return check_for_channel_wall_overlap();

  #endif

}

void collection_of_rigid_bodies::reset_to_start_of_step(){

  bodies = temp_bodies;

}

void collection_of_rigid_bodies::write_data_to_file(const double t) const {

  ofstream output_data_file(SIMULATION_NAME+std::string(".dat"),ios::app);
  output_data_file << std::scientific << std::setprecision(6);

  output_data_file << t << " ";

  for (int n=0; n<NUM_BODIES; n++){

    output_data_file << bodies[n].x(0) << " " << bodies[n].x(1) << " " << bodies[n].x(2) << " ";
    output_data_file << bodies[n].q.scalar_part << " " << bodies[n].q.vector_part(0) << " " << bodies[n].q.vector_part(1) << " " << bodies[n].q.vector_part(2) << " ";

  }

  output_data_file << endl;

  output_data_file.close();

  #if !MARKOV_CHAIN_MONTE_CARLO

  ofstream force_file(SIMULATION_NAME+std::string("_blob_forces.dat"),ios::app);
  force_file << std::scientific << std::setprecision(6);

  force_file << t << " ";

  for (int n=0; n<blob_forces.n_rows; n++){

    force_file << blob_forces(n) << " ";

  }

  force_file << endl;

  force_file.close();

  #endif

}

void collection_of_rigid_bodies::accept_configuration_from_master_process(const int myrank){

  double *temp = new double[7];

  for (int n=0; n<NUM_BODIES; n++){

    if (myrank==0){

      temp[0] = bodies[n].x(0);
      temp[1] = bodies[n].x(1);
      temp[2] = bodies[n].x(2);
      temp[3] = bodies[n].q.scalar_part;
      temp[4] = bodies[n].q.vector_part(0);
      temp[5] = bodies[n].q.vector_part(1);
      temp[6] = bodies[n].q.vector_part(2);

    }

    MPI_Bcast(temp, 7, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank!=0){

      bodies[n].accept_configuration_from_master_process(temp);

    }

  }

  delete[] temp;

}

double collection_of_rigid_bodies::calculate_divU(fcm_fluid_solver& fcm){

  double divU = 0.0;

  temp_bodies = bodies;

  #if SIMULATING_BLOBS

  for (int i=0; i<3; i++){

  #else

  for (int i=0; i<6; i++){

  #endif

    vec6 evec = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    if (i<3){

      evec(i) = FINITE_DIFFERENCE_DELTA_FACTOR*BLOB_RADIUS; // Scale by physical length scale for translational degrees of freedom.

    } else {

      evec(i) = FINITE_DIFFERENCE_DELTA_FACTOR;

    }

    for (int n=0; n<NUM_BODIES; n++){

      bodies[n].update(evec);

    }

    fill_all_blob_pos();
    fcm.set_up_gaussians(all_blob_pos);
    fcm.calculate_random_blob_velocities();

    vec V = fcm.V;

    for (int n=0; n<NUM_BODIES; n++){

      vec6 U0n = U0.col(n);

      #if SIMULATING_BLOBS

        divU += (V(3*n + i) - U0n(i))/evec(i);

      #else

        vec6 b = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

        const mat Q = bodies[n].q.matrix();

        for (int j=0; j<blobs_per_body; j++){

          int id = n*blobs_per_body + j;

          vec3 vel = V.subvec(3*id, 3*id + 2);

          b += join_vert(vel, cross(Q*blob_ref_pos.col(j), vel));

        }

        b.subvec(0,2) = Q.t() * b.subvec(0,2);
        b.subvec(3,5) = Q.t() * b.subvec(3,5);

        vec U = ref_KTKinv*b;

        U.subvec(0,2) = Q * U.subvec(0,2);
        U.subvec(3,5) = Q * U.subvec(3,5);

        divU += (U(i) - U0n(i))/evec(i);

      #endif

      bodies[n] = temp_bodies[n];

    }

  }

  return divU;

}
