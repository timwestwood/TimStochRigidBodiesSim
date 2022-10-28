#include <mpi.h>
#include <armadillo>
#include <iostream>
#include <iomanip>
#include "fcm_fluid_solver.hpp"
#include "collection_of_rigid_bodies.hpp"
#include "gmres_solver.hpp"
#include "config.hpp"

int main(int argc, char** argv){

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // 1. Start MPI
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  MPI_Init(&argc, &argv);

  int totalnodes, myrank;
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  for (int i=0; i<totalnodes; i++){

    if (i==myrank){

      cout << "Process " << i << "/" << totalnodes-1 << " is running..." << endl;

    }

    MPI_Barrier(MPI_COMM_WORLD);

  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // 2. Save parameters to file for post-processing
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (myrank == 0){

    ofstream par_file(SIMULATION_NAME+std::string(".par"));

    par_file << NUM_BODIES << " " << "%% NUM_BODIES" << endl;
    par_file << BLOB_RADIUS << " " << "%% BLOB_RADIUS" << endl;
    par_file << FLUID_VISCOSITY << " " << "%% FLUID_VISCOSITY" << endl;
    par_file << ENERGY_SCALE << " " << "%% ENERGY_SCALE" << endl;
    par_file << LX << " " << "%% LX" << endl;
    par_file << LY << " " << "%% LY" << endl;
    par_file << LZ << " " << "%% LZ" << endl;
    par_file << int(SLIP_CHANNEL)  << " " << "%% SLIP_CHANNEL" << endl;

    #if !MARKOV_CHAIN_MONTE_CARLO

    par_file << STEPS_PER_DIFFUSIVE_TIME << " " << "%% STEPS_PER_DIFFUSIVE_TIME" << endl;
    par_file << NUM_DIFFUSIVE_TIMES << " " << "%% NUM_DIFFUSIVE_TIMES" << endl;
    par_file << DATA_SAVES_PER_DIFFUSIVE_TIME << " " << "%% DATA_SAVES_PER_DIFFUSIVE_TIME" << endl;

    #if SHEARING_FORCES

    par_file << SHEAR_FORCE_SCALE << " " << "%% SHEAR_FORCE_SCALE" << endl;
    par_file << SHEAR_WAVENUMBER << " " << "%% SHEAR_WAVENUMBER" << endl;

    #endif

    #endif

    par_file.close();

  }

  MPI_Barrier(MPI_COMM_WORLD);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // 3. Prepare the solvers etc.
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  arma_rng::set_seed_random();

  collection_of_rigid_bodies all_bodies;
  all_bodies.initialise(myrank, totalnodes);

  #if SLIP_CHANNEL

  if (myrank == 0){
    cout << endl << "Volume fraction = " << 100.0*all_bodies.total_num_blobs*4.0*PI*BLOB_RADIUS*BLOB_RADIUS*BLOB_RADIUS/(3.0*LX*(0.5*LY)*LZ) << "%" << endl << endl;
  }

  #else

  if (myrank == 0){
    cout << endl << "Volume fraction = " << 100.0*all_bodies.total_num_blobs*4.0*PI*BLOB_RADIUS*BLOB_RADIUS*BLOB_RADIUS/(3.0*LX*LY*LZ) << "%" << endl << endl;
  }

  #endif

  #if !MARKOV_CHAIN_MONTE_CARLO

  fcm_fluid_solver fcm;
  fcm.initialise(all_bodies.total_num_blobs);
  fcm.produce_reference_matrices(all_bodies);

  #if SHEARING_FORCES

    fcm.produce_shear_velocity_field();

  #endif

  gmres_solver gmres;
  gmres.initialise(all_bodies.total_num_blobs);

  #endif

  all_bodies.seed_the_bodies(myrank);

  MPI_Barrier(MPI_COMM_WORLD);

  #if MARKOV_CHAIN_MONTE_CARLO

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // 4a. Begin sampling
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  all_bodies.fill_all_blob_pos();
  double U_old = all_bodies.calculate_potential();

  int n = 0;
  double scale = 0.1*BLOB_RADIUS;
  double acceptance_ratio = 0.5;
  const int ratio_steps = 20;
  while (n < 5000000){

    #if SIMULATING_BLOBS

    vec proposed_move(3*NUM_BODIES);

    if (myrank == 0){

      proposed_move = scale*(2.0*randu<vec>(3*NUM_BODIES) - 1.0);

    }

    MPI_Bcast(proposed_move.memptr(), 3*NUM_BODIES, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    #else

    vec proposed_move(6*NUM_BODIES);

    if (myrank == 0){

      proposed_move = scale*(2.0*randu<vec>(6*NUM_BODIES) - 1.0);

      for (int m=0; m<NUM_BODIES; m++){

        proposed_move.subvec(6*m + 3, 6*m + 5) /= all_bodies.body_radius;

      }

    }

    MPI_Bcast(proposed_move.memptr(), 6*NUM_BODIES, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    #endif

    if (all_bodies.end_of_step_update(proposed_move, 1.0)) {

      all_bodies.reset_to_start_of_step();
      continue;

    }

    double U = all_bodies.calculate_potential();

    if (randu<double>() < exp((U_old - U)/ENERGY_SCALE)){

      U_old = U;
      acceptance_ratio = ((ratio_steps - 1.0)*acceptance_ratio + 1.0)/double(ratio_steps);

    } else {

      all_bodies.reset_to_start_of_step();
      acceptance_ratio *= (ratio_steps - 1.0)/double(ratio_steps);

    }

    n++;

    if ((n%200 == 0) && (myrank == 0)){

      if (n > 100000){

        all_bodies.write_data_to_file(n);

        cout << "[scale/a = " << scale/BLOB_RADIUS << ", U = " << U << "] Saved sample " << n << endl;

      } else {

        cout << "[scale/a = " << scale/BLOB_RADIUS << ", U = " << U << "]" << endl;

      }

    }

    if (acceptance_ratio > 0.5){

      scale *= 1.02;

    } else {

      scale *= 0.98;

    }

    all_bodies.accept_configuration_from_master_process(myrank);

  }

  #else

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // 4b. Begin time-stepping
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  int skipped_steps = 0;
  int plotting_step = 0;
  double t = 0;
  const double final_time = NUM_DIFFUSIVE_TIMES * BLOB_DIFFUSIVE_TIME;

  while (t < final_time){

    if (myrank == 0){

      cout << "[t/tD = " << t/BLOB_DIFFUSIVE_TIME << "] Started the time-step..." << endl;

    }

    all_bodies.fill_all_blob_pos();
    fcm.produce_random_velocity_field();
    fcm.set_up_gaussians(all_bodies.all_blob_pos);
    fcm.calculate_random_blob_velocities();

    #if GEN_DRIFTER_CORRECTOR

      if (all_bodies.advance_to_midstep(fcm.V)) {

        if (myrank == 0){

          cout << "[t/tD = " << t/BLOB_DIFFUSIVE_TIME << "] Restarting the time-step due to a prohibited overlap..." << endl;

        }

        skipped_steps++;
        continue;

      }

      if (myrank == 0){

        cout << "[t/tD = " << t/BLOB_DIFFUSIVE_TIME << "] Reached the mid-step..." << endl;

      }

      const double divU = all_bodies.calculate_divU(fcm);

      const double nu = 1.0 + 0.5*DT*divU;

      if (myrank == 0){

        cout << "[t/tD = " << t/BLOB_DIFFUSIVE_TIME << "] Calculated div(U) = " << divU << endl;

      }

      all_bodies.fill_all_blob_pos_at_midstep();
      fcm.set_up_gaussians(all_bodies.all_blob_pos);
      fcm.calculate_random_blob_velocities();

      #if SHEARING_FORCES

        fcm.calculate_shear_blob_velocities();

      #endif

    #elif EULER_MARUYAMA

      const double nu = 1.0;

    #endif

    all_bodies.calculate_forces(t);

    if (myrank == 0){

      cout << "[t/tD = " << t/BLOB_DIFFUSIVE_TIME << "] Calculated the forces..." << endl;

    }

    #if SIMULATING_BLOBS

      vec V_stoch = fcm.V;

      #if SHEARING_FORCES

        V_stoch += fcm.Vshear;

      #endif

      fcm.set_up_force_distribution(all_bodies.blob_forces);
      fcm.produce_deterministic_velocity_field();
      fcm.calculate_deterministic_blob_velocities();

      bool overlap_at_end_of_step = all_bodies.end_of_step_update(V_stoch + fcm.V, nu);

    #else

      gmres.assemble_RHS(all_bodies, fcm);

      vec soln = gmres.invert_saddle_point_system(all_bodies, fcm, myrank, t);

      bool overlap_at_end_of_step = all_bodies.end_of_step_update(soln.subvec(3*all_bodies.total_num_blobs, 3*all_bodies.total_num_blobs + 6*NUM_BODIES - 1), nu);

    #endif

    if (overlap_at_end_of_step) {

      if (myrank == 0){

        cout << "[t/tD = " << t/BLOB_DIFFUSIVE_TIME << "] Restarting the time-step due to a prohibited overlap..." << endl;

      }

      all_bodies.reset_to_start_of_step();
      skipped_steps++;
      continue;

    }

    all_bodies.accept_configuration_from_master_process(myrank);

    t += DT;
    plotting_step++;

    if (plotting_step == ceil(STEPS_PER_DIFFUSIVE_TIME/DATA_SAVES_PER_DIFFUSIVE_TIME)){

      plotting_step = 0;

      if (myrank == 0){

        all_bodies.write_data_to_file(t);

        ofstream velocity_file(SIMULATION_NAME+std::string("_velocities.dat"),ios::app);
        velocity_file << std::scientific << std::setprecision(6);
        velocity_file << t << " ";

        #if SIMULATING_BLOBS

        for (int i=0; i<3*NUM_BODIES; i++){

          velocity_file << V_stoch(i) + fcm.V(i) << " ";

        }

        #else

        for (int i=0; i<6*NUM_BODIES; i++){

          velocity_file << soln(3*all_bodies.total_num_blobs + i) << " ";

        }

        #endif

        velocity_file << endl;

        velocity_file.close();

      }

      #if SHEARING_FORCES

      all_bodies.fill_all_blob_pos_at_midstep();
      fcm.set_up_gaussians(all_bodies.all_blob_pos);
      fcm.set_up_force_distribution(soln.subvec(0, 3*all_bodies.total_num_blobs - 1));
      fcm.produce_deterministic_velocity_field();

      ofstream fluid_velocity_file;

      if (myrank == 0){

        fluid_velocity_file.open(SIMULATION_NAME+std::string("_fluid_x_velocity_in_y.dat"), ios::app);

        fluid_velocity_file << t << " ";

      }

      const double *const ux = fcm.ux;
      const vec& ux_shear = fcm.ux_shear;

      for (int j = 0; j < NUM_FFT_GRID_POINTS_Y; j++){

        const int y_index = j*fcm.fft.pad;

        double u_local = 0.0;

        for (int i = 0; i < fcm.fft.local_nx; i++){

          const int x_index = i*NUM_FFT_GRID_POINTS_Y*fcm.fft.pad;

          for (int k = 0; k < NUM_FFT_GRID_POINTS_Z; k++){

            const int index = x_index + y_index + k;

            u_local += ux[index] + ux_shear(index); // ux_rand should average to zero over time

          }

        }

        double u_global = 0.0;

        MPI_Reduce(&u_local, &u_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); // Sum onto process 0 to write to file.

        if (myrank == 0){

          u_global /= NUM_FFT_GRID_POINTS_X*NUM_FFT_GRID_POINTS_Z; // Average across the entire x and z domains.

          fluid_velocity_file << u_global << " "; // j^th entry on a given row of the file corresponds to y = (j-1)*DX

        }

      }

      if (myrank == 0){

        fluid_velocity_file << endl;

        fluid_velocity_file.close();

      }

      #endif

      if (myrank == 0){

        cout << "Reached time t/tD = " << t/BLOB_DIFFUSIVE_TIME << " and saved to file." << endl;

      }

    } else {

      if (myrank == 0){

        cout << "Reached time t/tD = " << t/BLOB_DIFFUSIVE_TIME << "." << endl;

      }

    }

  }

  if (myrank == 0){

    ofstream par_file(SIMULATION_NAME+std::string(".par"),ios::app);
    par_file << skipped_steps << " " << "%% Number of steps skipped due to prohibited overlapping." << endl;
    par_file.close();

  }

  #endif

  MPI_Finalize();

  return 0;

}
