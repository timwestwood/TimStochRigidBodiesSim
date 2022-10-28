// gmres_solver.cpp

#include <iostream>
#include "collection_of_rigid_bodies.hpp"
#include "fcm_fluid_solver.hpp"
#include "gmres_solver.hpp"
#include "rcross.hpp"
#include "config.hpp"

gmres_solver::~gmres_solver(){}

gmres_solver::gmres_solver(){}

void gmres_solver::initialise(const int total_num_blobs){

  num_blobs = total_num_blobs;

  system_size = 3*num_blobs + 6*NUM_BODIES;

  RHS = vec(system_size,fill::zeros);

  const int max_iter = std::min<int>(system_size,GMRES_MAX_ITER);

  Q = mat(system_size, max_iter+1, fill::zeros);

  beta = vec(max_iter+1, fill::zeros);

  H = mat(max_iter+1, max_iter+1, fill::zeros);

  SN = vec(max_iter, fill::zeros);

  CS = vec(max_iter, fill::zeros);

}

/*

The system we want to solve is

      K*V = M*lambda + Ub
      K^T * lambda = -F

where M is the mobility matrix for all of the blobs comprising rigid bodies, Ub contains the random velocities
of all of the blobs, F contains the total body forces and torques, K is the Jacobian of the rigid
body constraints, lambda contains the forces on the blobs (including the Lagrange multipliers associated
with the rigid body constraints) and V contains the translational and angular velocities of the rigid bodies.

Written in matrix form, this system of equations becomes

    |  M     -K |     | lambda |      | -Ub |
    |           |  *  |        |  =   |     |
    |-K^T     0 |     |   V    |      | -F  |

The purpose of this GMRES class is to solve this system for V (and lambda, but we don't actually use it).

*/

void gmres_solver::assemble_RHS(const collection_of_rigid_bodies& all_bodies, const fcm_fluid_solver& fcm){

  RHS.subvec(0, 3*num_blobs-1) = -fcm.V;

  #if SHEARING_FORCES

    RHS.subvec(0, 3*num_blobs-1) -= fcm.Vshear;

  #endif

  RHS.subvec(3*num_blobs, system_size-1) = -all_bodies.multiply_by_KT_at_midstep(all_bodies.blob_forces);

}

vec gmres_solver::multiply_by_saddle_point_matrix(const vec& in, const collection_of_rigid_bodies& all_bodies, fcm_fluid_solver& fcm){

  fcm.set_up_force_distribution(in.subvec(0,3*num_blobs-1));
  fcm.produce_deterministic_velocity_field();
  fcm.calculate_deterministic_blob_velocities();
  vec out = fcm.V;

  out -= all_bodies.multiply_by_K_at_midstep(in.subvec(3*num_blobs,system_size-1));

  return join_vert(out,-all_bodies.multiply_by_KT_at_midstep(in.subvec(0,3*num_blobs-1)));

}

vec gmres_solver::multiply_by_Pinv(const vec& in, const collection_of_rigid_bodies& all_bodies, const fcm_fluid_solver& fcm){

  /*

  This function multiplies a given vector by the matrix P^(-1) where P is an approximation to the saddle-point
  system matrix. To do this, we make a block-diagonal approximation to the blob mobility matrix and then invert the
  approximate system using Schur complements.

  */

  vec out(system_size);

  mat K(3*all_bodies.blobs_per_body, 6);

  for (int m=0; m < all_bodies.blobs_per_body; m++){

    K(3*m, 0, size(3,3)) = eye(3,3);

  }

  for (int n=0; n<NUM_BODIES; n++){

    const mat Q = all_bodies.midstep_bodies[n].q.matrix();

    // Form K

    for (int m=0; m < all_bodies.blobs_per_body; m++){

      K(3*m, 3, size(3,3)) = rcross(Q * all_bodies.blob_ref_pos.col(m));

    }

    // Map the reference matrices to the current configuration

    mat Minv = fcm.ref_Minv;
    mat N = fcm.ref_N;

    for (int m=0; m < all_bodies.blobs_per_body; m++){
      for (int l=m; l < all_bodies.blobs_per_body; l++){

        Minv(3*m, 3*l, size(3,3)) = Q*Minv(3*m, 3*l, size(3,3))*Q.t();
        Minv(3*l, 3*m, size(3,3)) = Minv(3*m, 3*l, size(3,3)).t();

      }
    }

    N(0, 0, size(3,3)) = Q*N(0, 0, size(3,3))*Q.t();
    N(0, 3, size(3,3)) = Q*N(0, 3, size(3,3))*Q.t();
    N(3, 0, size(3,3)) = Q*N(3, 0, size(3,3))*Q.t();
    N(3, 3, size(3,3)) = Q*N(3, 3, size(3,3))*Q.t();

    // Solve the block of the system for body n
    out.subvec(3*num_blobs + 6*n, 3*num_blobs + 6*n + 5) = -in.subvec(3*num_blobs + 6*n, 3*num_blobs + 6*n + 5) - K.t()*Minv*in.subvec(3*all_bodies.blobs_per_body*n, 3*all_bodies.blobs_per_body*(n+1) - 1);
    out.subvec(3*num_blobs + 6*n, 3*num_blobs + 6*n + 5) = N * out.subvec(3*num_blobs + 6*n, 3*num_blobs + 6*n + 5);

    out.subvec(3*all_bodies.blobs_per_body*n, 3*all_bodies.blobs_per_body*(n+1) - 1) = Minv*(K*out.subvec(3*num_blobs + 6*n, 3*num_blobs + 6*n + 5) + in.subvec(3*all_bodies.blobs_per_body*n, 3*all_bodies.blobs_per_body*(n+1) - 1));

  }

  return out;

}

vec gmres_solver::invert_saddle_point_system(const collection_of_rigid_bodies& all_bodies, fcm_fluid_solver& fcm, const int myrank, const double t){

  /*

  In this function, we solve the system A * x = RHS, where A is the saddle-point matrix or a preconditioned version of it.

  We search iteratively for a solution x_n in the n-th Krylov subspace; i.e. the space spanned by RHS, A * RHS, ..., A^(n-1) * RHS.

  Since these spanning vectors may be close to linearly dependent, we use the stabilised Gram-Schmidt process to construct an orthonormal basis
  q_1, ..., q_n, which we store as the columns of the (system_size)-by-n matrix Q_n. Hence, there exists some vector y_n of coefficients such that
  x_n = Q_n * y_n.

  As we produce this orthonormal basis, we can also construct the (n+1)-by-n upper Hessenburg matrix H_n which satisfies A * Q_n = Q_(n+1) * H_n;
  the first n elements in column n are the projections of the corresponding Gram-Schmidt iterate onto q_n, and the (n+1)-th element is the norm of
  q_(n+1) before it is projected to unit length. The previous columns are unchanged -- we form H_n from H_(n-1) by adding a row of zeros and the
  column just described.

  Since the columns of Q_n are orthonormal,

    error = ||A * x_n - RHS|| = ||H_n * y_n - Q^T_(n+1) * RHS || = ||H_n * y_n - ||RHS|| * e_1 ||.

  Thus we can find x_n by solving a least-squares problem for y_n. Although H_n is not a square matrix, we can convert the problem to a square one
  as we proceed by performing successive Givens rotations to both the H matrices and the vector which starts as ||RHS|| * e_1. This process will
  also give us the error in the least-squares problem, and hence the error in the overall linear system.

  N.B. We actually solve the equivalent pre-conditioned system A*P^(-1)*P*x = RHS, where P is a readily invertible approximation to A.

  */

  Q.zeros();

  const double norm_of_RHS = norm(RHS);

  Q.col(0) = RHS/norm_of_RHS;

  beta.zeros();

  beta(0) = norm_of_RHS;

  H.zeros();
  CS.zeros();
  SN.zeros();

  const int max_iter = std::min<int>(system_size,GMRES_MAX_ITER);

  for (int iter = 1; iter <= max_iter; iter++){

    // Produce the new orthonormal vector, using the appropriate values to update H as we do.
    Q.col(iter) = multiply_by_saddle_point_matrix(multiply_by_Pinv(Q.col(iter-1), all_bodies, fcm), all_bodies, fcm);

    for (int i=0; i<iter; i++){

      H(i,iter-1) = dot(Q.col(iter),Q.col(i));

      Q.col(iter) -= H(i,iter-1)*Q.col(i);

    }

    double q_norm = norm(Q.col(iter));

    H(iter,iter-1) = q_norm;

    Q.col(iter) /= q_norm;

    // Apply Givens rotations to transform the least-squares problem to a square one.
    // Apply any previous rotations IN ORDER. These have all already been applied to the RHS vector beta.
    for (int i=1; i<iter; i++){

      double temp = CS(i-1)*H(i-1,iter-1) + SN(i-1)*H(i,iter-1);
      H(i,iter-1) = CS(i-1)*H(i,iter-1) - SN(i-1)*H(i-1,iter-1);
      H(i-1,iter-1) = temp;

    }

    // Then calculate and apply the new Givens rotation. This one is also applied to the RHS of the least-squares problem, beta.

    if (H(iter,iter-1)==0.0){

      // It's already lower diagonal.
      CS(iter-1) = 1.0;
      SN(iter-1) = 0.0;

    } else {

      if (fabs(H(iter,iter-1)) > fabs(H(iter-1,iter-1))){

        double temp = H(iter-1,iter-1)/H(iter,iter-1);
        SN(iter-1) = (2.0*double(H(iter,iter-1) > 0) - 1.0)/sqrt(1.0 + temp*temp);
        CS(iter-1) = temp*SN(iter-1);

      } else {

        double temp = H(iter,iter-1)/H(iter-1,iter-1);
        CS(iter-1) = (2.0*double(H(iter-1,iter-1) > 0) - 1.0)/sqrt(1.0 + temp*temp);
        SN(iter-1) = temp*CS(iter-1);

      }

    }

    H(iter-1,iter-1) = CS(iter-1)*H(iter-1,iter-1) + SN(iter-1)*H(iter,iter-1);
    H(iter,iter-1) = 0.0;

    beta(iter) = -SN(iter-1)*beta(iter-1);
    beta(iter-1) = CS(iter-1)*beta(iter-1);

    // beta(iter) now contains the (signed) error in the least-squares system, and hence the error in the linear system.
    // If it is small enough, or we have hit the maximum number of iterations, we generate the solution and return.
    const double relative_error = fabs(beta(iter))/norm_of_RHS;

    if ((relative_error < GMRES_TOL) || (iter == max_iter)){

      vec y = vec(iter,fill::zeros);

      y(iter-1) = beta(iter-1)/H(iter-1,iter-1);

      for (int i=iter-2; i>=0; i--){

        y(i) = beta(i);

        for (int j=iter-1; j>i; j--){

          y(i) -= H(i,j)*y(j);

        }

        y(i) /= H(i,i);

      }

      vec soln = vec(system_size,fill::zeros);

      for (int i=0; i<iter; i++){

        soln += y(i)*Q.col(i);

      }

      if (myrank == 0){

        if (relative_error < GMRES_TOL){

          cout << "[t/tD = " << t/BLOB_DIFFUSIVE_TIME << "] GMRES converged to a relative error of " << relative_error << " in " << iter << " iterations..." << endl;

        } else {

          cout << "[t/tD = " << t/BLOB_DIFFUSIVE_TIME << "] GMRES failed to converge in " << iter << " iterations. The relative error is " << relative_error << endl;

        }

      }

      return multiply_by_Pinv(soln, all_bodies, fcm);

    }

  }

}
