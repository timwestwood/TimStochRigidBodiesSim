// fcm_fluid_solver.cpp

#include <mpi.h>
#include "parallel_random_numbers.hpp"
#include "collection_of_rigid_bodies.hpp"
#include "fcm_fluid_solver.hpp"
#include "rcross.hpp"
#include "config.hpp"

fcm_fluid_solver::~fcm_fluid_solver(){

  delete[] ux;
  delete[] uy;
  delete[] uz;
  delete[] fx;
  delete[] fy;
  delete[] fz;
  delete[] sxx;
  delete[] syy;
  delete[] szz;
  delete[] sxy;
  delete[] sxz;
  delete[] syz;

}

fcm_fluid_solver::fcm_fluid_solver(){}

void fcm_fluid_solver::initialise(const int total_num_blobs){

  fft.initialise();

  ux = new double[fft.total_local_size];
  uy = new double[fft.total_local_size];
  uz = new double[fft.total_local_size];
  fx = new double[fft.total_local_size];
  fy = new double[fft.total_local_size];
  fz = new double[fft.total_local_size];
  sxx = new double[fft.total_local_size];
  syy = new double[fft.total_local_size];
  szz = new double[fft.total_local_size];
  sxy = new double[fft.total_local_size];
  sxz = new double[fft.total_local_size];
  syz = new double[fft.total_local_size];

  V = vec(3*total_num_blobs);
  Vtemp = vec(3*total_num_blobs);

  ux_rand = vec(fft.total_local_size);
  uy_rand = vec(fft.total_local_size);
  uz_rand = vec(fft.total_local_size);

  #if SHEARING_FORCES

    Vshear = vec(3*total_num_blobs);

    ux_shear = vec(fft.total_local_size);
    uy_shear = vec(fft.total_local_size);
    uz_shear = vec(fft.total_local_size);

  #endif

  gaussx = mat(total_num_blobs, NUM_FCM_GAUSSIAN_POINTS);
  gaussy = mat(total_num_blobs, NUM_FCM_GAUSSIAN_POINTS);
  gaussz = mat(total_num_blobs, NUM_FCM_GAUSSIAN_POINTS);

  indx = Mat<int>(total_num_blobs, NUM_FCM_GAUSSIAN_POINTS);
  indy = Mat<int>(total_num_blobs, NUM_FCM_GAUSSIAN_POINTS);
  indz = Mat<int>(total_num_blobs, NUM_FCM_GAUSSIAN_POINTS);

  qx = vec(NUM_FFT_GRID_POINTS_X);
  qxsq = vec(NUM_FFT_GRID_POINTS_X);

  int nptsh =  (NUM_FFT_GRID_POINTS_X / 2);
  double fourier_fac = 2.0*PI/LX;

   for (int i = 0; i < NUM_FFT_GRID_POINTS_X; i++){

    if (i <= nptsh){

      qx(i) = i*fourier_fac;

    } else {

      qx(i) = (i - NUM_FFT_GRID_POINTS_X)*fourier_fac;

    }

    qxsq(i) = qx(i)*qx(i);

  }

  qy = vec(NUM_FFT_GRID_POINTS_Y);
  qysq = vec(NUM_FFT_GRID_POINTS_Y);

  nptsh =  (NUM_FFT_GRID_POINTS_Y / 2);
  fourier_fac = 2.0*PI/LY;

   for (int i = 0; i < NUM_FFT_GRID_POINTS_Y; i++){

    if (i <= nptsh){

      qy(i) = i*fourier_fac;

    } else {

      qy(i) = (i - NUM_FFT_GRID_POINTS_Y)*fourier_fac;

    }

    qysq(i) = qy(i)*qy(i);

  }

  qpad = vec(fft.pad);
  qpadsq = vec(fft.pad);

  fourier_fac = 2.0*PI/LZ;

  for (int i = 0; i < fft.pad; i = i+2){

    qpad(i) = 0.5*i*fourier_fac;
    qpad(i+1) = qpad(i);
    qpadsq(i) = qpad(i) * qpad(i);
    qpadsq(i+1) = qpadsq(i);

  }

}

void fcm_fluid_solver::set_up_gaussians(const mat& Y){

  const int ngdh = NUM_FCM_GAUSSIAN_POINTS/2;

  const double anorm = 1.0/(sqrt(2.0)*BLOB_RADIUS); // See eqns (4) and (5) of Lomholt and Maxey (2003).
  const double anorm2 = 2.0*BLOB_RADIUS*BLOB_RADIUS/PI;

  for (int np = 0; np < Y.n_rows; np++){

    const int xc = int(Y(np,0)/DX);
    const int yc = int(Y(np,1)/DX);
    const int zc = int(Y(np,2)/DX);

    for (int i = 0; i < NUM_FCM_GAUSSIAN_POINTS; i++){

      int xg = xc - ngdh + (i+1);
      indx(np,i) = xg - NUM_FFT_GRID_POINTS_X * int(floor(double(xg)/double(NUM_FFT_GRID_POINTS_X)));
      double xx = xg*DX - Y(np,0);
      gaussx(np,i) = anorm*exp(-xx*xx/anorm2);

      int yg = yc - ngdh + (i+1);
      indy(np,i) = yg - NUM_FFT_GRID_POINTS_Y * int(floor(double(yg)/double(NUM_FFT_GRID_POINTS_Y)));
      xx = yg*DX - Y(np,1);
      gaussy(np,i) = anorm*exp(-xx*xx/anorm2);

      int zg = zc - ngdh + (i+1);
      indz(np,i) = zg - NUM_FFT_GRID_POINTS_Z * int(floor(double(zg)/double(NUM_FFT_GRID_POINTS_Z)));
      xx = zg*DX - Y(np,2);
      gaussz(np,i) = anorm*exp(-xx*xx/anorm2);

    }

  }

}

void fcm_fluid_solver::produce_random_velocity_field(){

  const double var1 = sqrt(2.0*ENERGY_SCALE*FLUID_VISCOSITY/(DX*DX*DX*DT));
  const double var2 = sqrt(2.0)*var1;

  const int half_y = NUM_FFT_GRID_POINTS_Y/2;

  for (int i = 0; i < fft.local_nx; i++){

    const int x_index = i*NUM_FFT_GRID_POINTS_Y*fft.pad;

    for (int k = 0; k < fft.pad; k++){

      const int z_index = k;

      for (int j=0; j<NUM_FFT_GRID_POINTS_Y; j++){

        const int index = x_index + z_index + j*fft.pad;

        #if SLIP_CHANNEL

          if ((j == 0) || (j == half_y)){

            sxx[index] = sqrt(2.0)*var2*gaussian_random();
            syy[index] = sqrt(2.0)*var2*gaussian_random();
            szz[index] = sqrt(2.0)*var2*gaussian_random();
            sxz[index] = sqrt(2.0)*var1*gaussian_random();
            sxy[index] = 0.0;
            syz[index] = 0.0;

          } else if (j < half_y){

            sxx[index] = var2*gaussian_random();
            syy[index] = var2*gaussian_random();
            szz[index] = var2*gaussian_random();
            sxz[index] = var1*gaussian_random();
            sxy[index] = var1*gaussian_random();
            syz[index] = var1*gaussian_random();

          } else {

            const int image_index = x_index + z_index + (NUM_FFT_GRID_POINTS_Y-j)*fft.pad;

            sxx[index] = sxx[image_index];
            syy[index] = syy[image_index];
            szz[index] = szz[image_index];
            sxz[index] = sxz[image_index];
            sxy[index] = -sxy[image_index];
            syz[index] = -syz[image_index];

          }

        #else

          sxx[index] = var2*gaussian_random();
          syy[index] = var2*gaussian_random();
          szz[index] = var2*gaussian_random();
          sxy[index] = var1*gaussian_random();
          sxz[index] = var1*gaussian_random();
          syz[index] = var1*gaussian_random();

        #endif

      }

    }

  }

  fft.transform(sxx);
  fft.transform(sxy);
  fft.transform(sxz);
  fft.transform(syy);
  fft.transform(szz);
  fft.transform(syz);

  for (int j = 0; j < fft.local_ny_after_transpose; j++){

    double q2 = qy(j + fft.local_y_start_after_transpose);

    if (j + fft.local_y_start_after_transpose == NUM_FFT_GRID_POINTS_Y/2){

      q2 = 0.0;

    }

    const int y_index = j*NUM_FFT_GRID_POINTS_X*fft.pad;

    for (int i = 0; i < NUM_FFT_GRID_POINTS_X; i++){

      double q1 = qx(i);

      if (i == NUM_FFT_GRID_POINTS_X/2){

        q1 = 0.0;

      }

      const int x_index = i*fft.pad;

      for (int k = 0; k < fft.pad; k=k+2){

        const int index = x_index + y_index + k;
        const int index_plus_one = index + 1;

        double q3 = qpad(k);

        if (k == fft.pad-2){

          q3 = 0.0;

        }

        // e.g. fx[index] + sqrt(-1)*fx[index_plus_one] = fx(q1,q2,q3) = sqrt(-1)*(q1*sxx(q1,q2,q3) + q2*sxy(q1,q2,q3) + q3*sxz(q1,q2,q3)) = sqrt(-1)*(q1*(sxx[index] + sqrt(-1)*sxx[index_plus_one]) + ...)

        fx[index] = -q1*sxx[index_plus_one] - q2*sxy[index_plus_one] - q3*sxz[index_plus_one];
        fy[index] = -q1*sxy[index_plus_one] - q2*syy[index_plus_one] - q3*syz[index_plus_one];
        fz[index] = -q1*sxz[index_plus_one] - q2*syz[index_plus_one] - q3*szz[index_plus_one];

        fx[index_plus_one] = q1*sxx[index] + q2*sxy[index] + q3*sxz[index];
        fy[index_plus_one] = q1*sxy[index] + q2*syy[index] + q3*syz[index];
        fz[index_plus_one] = q1*sxz[index] + q2*syz[index] + q3*szz[index];

      }
    }
  }

  const double totalpts = double(NUM_FFT_GRID_POINTS_X * NUM_FFT_GRID_POINTS_Y * NUM_FFT_GRID_POINTS_Z);

  if (fft.local_y_start_after_transpose == 0){

    fx[0] = 0.0;
    fx[1] = 0.0;
    fy[0] = 0.0;
    fy[1] = 0.0;
    fz[0] = 0.0;
    fz[1] = 0.0;

  }

  for (int j = 0; j < fft.local_ny_after_transpose; j++){

    double q2 = qy(j + fft.local_y_start_after_transpose);

    const int y_index = j*NUM_FFT_GRID_POINTS_X*fft.pad;

    for (int i = 0; i < NUM_FFT_GRID_POINTS_X; i++){

      double q1 = qx(i);
      double qq = qysq(j + fft.local_y_start_after_transpose) + qxsq(i) + 1e-16;

      const int x_index = i*fft.pad;

      for (int k = 0; k < fft.pad; k++){

        const int index = x_index + y_index + k;

        const double q3 = qpad(k);
        const double f1 = fx[index];
        const double f2 = fy[index];
        const double f3 = fz[index];
        const double norm = 1.0/(qq + qpadsq(k));
        const double kdotf = (q1*f1+q2*f2+q3*f3)*norm;
        ux[index] = norm*(f1-q1*(kdotf))/(FLUID_VISCOSITY*totalpts);
        uy[index] = norm*(f2-q2*(kdotf))/(FLUID_VISCOSITY*totalpts);
        uz[index] = norm*(f3-q3*(kdotf))/(FLUID_VISCOSITY*totalpts);
        
      }
    }
  }

  if (fft.local_y_start_after_transpose == 0){

    ux[0] = 0.0;
    ux[1] = 0.0;
    uy[0] = 0.0;
    uy[1] = 0.0;
    uz[0] = 0.0;
    uz[1] = 0.0;

  }

  fft.inverse_transform(ux);
  fft.inverse_transform(uy);
  fft.inverse_transform(uz);

  for (int index = 0; index < fft.total_local_size; index++){

    ux_rand(index) = ux[index];
    uy_rand(index) = uy[index];
    uz_rand(index) = uz[index];

  }

}

void fcm_fluid_solver::calculate_deterministic_blob_velocities(){

  const int num_blobs = gaussx.n_rows;

  Vtemp.zeros();

  const double norm = DX*DX*DX;

  for (int np = 0; np < num_blobs; np++){

    for (int i = 0; i < NUM_FCM_GAUSSIAN_POINTS; i++){

      int ii = indx(np,i);

      if ((fft.local_x_start <= ii) && (ii < (fft.local_x_start + fft.local_nx))){

         ii -= fft.local_x_start;

         const int x_index = ii*NUM_FFT_GRID_POINTS_Y*fft.pad;

         for (int j = 0; j < NUM_FCM_GAUSSIAN_POINTS; j++){

            const int jj = indy(np,j);

            #if SLIP_CHANNEL

            if (jj < NUM_FFT_GRID_POINTS_Y/2){ // Truncate Gaussians

            #endif

            const int y_index = jj*fft.pad;

            for (int k = 0; k < NUM_FCM_GAUSSIAN_POINTS; k++){

               const int kk = indz(np,k);

               const int index = x_index + y_index + kk;

               const double temp = gaussx(np,i)*gaussy(np,j)*gaussz(np,k)*norm;

               Vtemp(3*np) += ux[index]*temp;
               Vtemp(3*np + 1) += uy[index]*temp;
               Vtemp(3*np + 2) += uz[index]*temp;

              }

            #if SLIP_CHANNEL

            }

            #endif

             }
           }
         }
       }

  V.zeros();

  MPI_Allreduce(Vtemp.memptr(), V.memptr(), 3*num_blobs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}

void fcm_fluid_solver::calculate_random_blob_velocities(){

  const int num_blobs = gaussx.n_rows;

  Vtemp.zeros();

  const double norm = DX*DX*DX;

  for (int np = 0; np < num_blobs; np++){

    for (int i = 0; i < NUM_FCM_GAUSSIAN_POINTS; i++){

      int ii = indx(np,i);

      if ((fft.local_x_start <= ii) && (ii < (fft.local_x_start + fft.local_nx))){

         ii -= fft.local_x_start;

         const int x_index = ii*NUM_FFT_GRID_POINTS_Y*fft.pad;

         for (int j = 0; j < NUM_FCM_GAUSSIAN_POINTS; j++){

            const int jj = indy(np,j);

            #if SLIP_CHANNEL

            if (jj < NUM_FFT_GRID_POINTS_Y/2){ // Truncate Gaussians

            #endif

            const int y_index = jj*fft.pad;

            for (int k = 0; k < NUM_FCM_GAUSSIAN_POINTS; k++){

               const int kk = indz(np,k);

               const int index = x_index + y_index + kk;

               const double temp = gaussx(np,i)*gaussy(np,j)*gaussz(np,k)*norm;

               Vtemp(3*np) += ux_rand(index)*temp;
               Vtemp(3*np + 1) += uy_rand(index)*temp;
               Vtemp(3*np + 2) += uz_rand(index)*temp;

              }

            #if SLIP_CHANNEL

            }

            #endif

             }
           }
         }
       }

  V.zeros();

  MPI_Allreduce(Vtemp.memptr(), V.memptr(), 3*num_blobs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}

void fcm_fluid_solver::calculate_shear_blob_velocities(){

  const int num_blobs = gaussx.n_rows;

  Vtemp.zeros();

  const double norm = DX*DX*DX;

  for (int np = 0; np < num_blobs; np++){

    for (int i = 0; i < NUM_FCM_GAUSSIAN_POINTS; i++){

      int ii = indx(np,i);

      if ((fft.local_x_start <= ii) && (ii < (fft.local_x_start + fft.local_nx))){

         ii -= fft.local_x_start;

         const int x_index = ii*NUM_FFT_GRID_POINTS_Y*fft.pad;

         for (int j = 0; j < NUM_FCM_GAUSSIAN_POINTS; j++){

            const int jj = indy(np,j);

            #if SLIP_CHANNEL

            if (jj < NUM_FFT_GRID_POINTS_Y/2){ // Truncate Gaussians

            #endif

            const int y_index = jj*fft.pad;

            for (int k = 0; k < NUM_FCM_GAUSSIAN_POINTS; k++){

               const int kk = indz(np,k);

               const int index = x_index + y_index + kk;

               const double temp = gaussx(np,i)*gaussy(np,j)*gaussz(np,k)*norm;

               Vtemp(3*np) += ux_shear(index)*temp;
               Vtemp(3*np + 1) += uy_shear(index)*temp;
               Vtemp(3*np + 2) += uz_shear(index)*temp;

              }

            #if SLIP_CHANNEL

            }

            #endif

             }
           }
         }
       }

  Vshear.zeros();

  MPI_Allreduce(Vtemp.memptr(), Vshear.memptr(), 3*num_blobs, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

}

void fcm_fluid_solver::set_up_force_distribution(const vec& F){

  const int Nblobs = gaussx.n_rows;

  for (int index = 0; index < fft.total_local_size; index++){

    fx[index] = 0.0;
    fy[index] = 0.0;
    fz[index] = 0.0;

  }

  for (int np = 0; np < Nblobs; np++){

    for (int i = 0; i < NUM_FCM_GAUSSIAN_POINTS; i++){

      int ii = indx(np,i);

      if ((fft.local_x_start <= ii) && (ii < (fft.local_x_start + fft.local_nx))){

	      ii -= fft.local_x_start;

        const int x_index = ii*NUM_FFT_GRID_POINTS_Y*fft.pad;

        for (int j = 0; j < NUM_FCM_GAUSSIAN_POINTS; j++){

          int jj = indy(np,j);

          const int y_index = jj*fft.pad;

          for (int k = 0; k < NUM_FCM_GAUSSIAN_POINTS; k++){

             int kk = indz(np,k);

             const int index = x_index + y_index + kk;

             double temp = gaussx(np,i)*gaussy(np,j)*gaussz(np,k);

             int id = 3*np;

        	   fx[index] += F(id)*temp;
        	   fy[index] += F(id+1)*temp;
        	   fz[index] += F(id+2)*temp;

          }
        }
      }
    }
  }

  #if SLIP_CHANNEL

  const int half_y = NUM_FFT_GRID_POINTS_Y/2;

  for (int i=0; i < fft.local_nx; i++){

    const int x_index = i*NUM_FFT_GRID_POINTS_Y*fft.pad;

    for (int k=0; k < fft.pad; k++){

      const int z_index = k;

      fy[x_index + z_index] = 0.0; // fy[i][0][k] = 0.0;
      fy[x_index + z_index + half_y*fft.pad] = 0.0; // fy[i][half_y][k] = 0.0;

      for (int j = (1 + half_y); j < NUM_FFT_GRID_POINTS_Y; j++){

        const int index = x_index + z_index + j*fft.pad;
        const int image_index = x_index + z_index + (NUM_FFT_GRID_POINTS_Y-j)*fft.pad;

        fx[index] = fx[image_index];
        fy[index] = -fy[image_index];
        fz[index] = fz[image_index];

      }

    }

  }

  #endif

}

void fcm_fluid_solver::produce_deterministic_velocity_field(){

  fft.transform(fx);
  fft.transform(fy);
  fft.transform(fz);

  if (fft.local_y_start_after_transpose == 0){

    fx[0] = 0.0;
    fx[1] = 0.0;
    fy[0] = 0.0;
    fy[1] = 0.0;
    fz[0] = 0.0;
    fz[1] = 0.0;

  }

  const double totalpts = double(NUM_FFT_GRID_POINTS_X * NUM_FFT_GRID_POINTS_Y * NUM_FFT_GRID_POINTS_Z);

  for (int j = 0; j < fft.local_ny_after_transpose; j++){

    const double q2 = qy(j + fft.local_y_start_after_transpose);

    const int y_index = j*NUM_FFT_GRID_POINTS_X*fft.pad;

    for (int i = 0; i < NUM_FFT_GRID_POINTS_X; i++){

      const double q1 = qx(i);
      const double qq = qysq(j + fft.local_y_start_after_transpose) + qxsq(i) + 1e-16;

      const int x_index = i*fft.pad;

      for (int k = 0; k < fft.pad; k++){

        const int index = x_index + y_index + k;

        const double q3 = qpad(k);
        const double f1 = fx[index];
        const double f2 = fy[index];
        const double f3 = fz[index];
        const double norm = 1.0/(qq + qpadsq(k));
        const double kdotf = (q1*f1+q2*f2+q3*f3)*norm;
        ux[index] = norm*(f1-q1*(kdotf))/(FLUID_VISCOSITY*totalpts);
        uy[index] = norm*(f2-q2*(kdotf))/(FLUID_VISCOSITY*totalpts);
        uz[index] = norm*(f3-q3*(kdotf))/(FLUID_VISCOSITY*totalpts);

      }
    }
  }

  if (fft.local_y_start_after_transpose == 0){

    ux[0] = 0.0;
    ux[1] = 0.0;
    uy[0] = 0.0;
    uy[1] = 0.0;
    uz[0] = 0.0;
    uz[1] = 0.0;

  }

  fft.inverse_transform(ux);
  fft.inverse_transform(uy);
  fft.inverse_transform(uz);

}

void fcm_fluid_solver::produce_shear_velocity_field(){

  #if SHEARING_FORCES

  for (int index = 0; index < fft.total_local_size; index++){

    fx[index] = 0.0;
    fy[index] = 0.0;
    fz[index] = 0.0;

  }

  for (int i = 0; i < fft.local_nx; i++){

    const int x_index = i*NUM_FFT_GRID_POINTS_Y*fft.pad;

    for (int j = 0; j < NUM_FFT_GRID_POINTS_Y; j++){

      const int y_index = j*fft.pad;

      for (int k = 0; k < NUM_FFT_GRID_POINTS_Z; k++){

        const int index = x_index + y_index + k;

        fx[index] = SHEAR_FORCE_SCALE*ENERGY_SCALE*sin(SHEAR_WAVENUMBER*2.0*PI*double(j)/double(NUM_FFT_GRID_POINTS_Y))/BLOB_RADIUS;

      }

    }

  }

  #if SLIP_CHANNEL

  const int half_y = NUM_FFT_GRID_POINTS_Y/2;

  for (int i=0; i < fft.local_nx; i++){

    const int x_index = i*NUM_FFT_GRID_POINTS_Y*fft.pad;

    for (int k=0; k < fft.pad; k++){

      const int z_index = k;

      fy[x_index + z_index] = 0.0; // fy[i][0][k] = 0.0;
      fy[x_index + z_index + half_y*fft.pad] = 0.0; // fy[i][half_y][k] = 0.0;

      for (int j = (1 + half_y); j < NUM_FFT_GRID_POINTS_Y; j++){

        const int index = x_index + z_index + j*fft.pad;
        const int image_index = x_index + z_index + (NUM_FFT_GRID_POINTS_Y-j)*fft.pad;

        fx[index] = fx[image_index];
        fy[index] = -fy[image_index];
        fz[index] = fz[image_index];

      }

    }

  }

  #endif

  fft.transform(fx);
  fft.transform(fy);
  fft.transform(fz);

  if (fft.local_y_start_after_transpose == 0){

    fx[0] = 0.0;
    fx[1] = 0.0;
    fy[0] = 0.0;
    fy[1] = 0.0;
    fz[0] = 0.0;
    fz[1] = 0.0;

  }

  const double totalpts = double(NUM_FFT_GRID_POINTS_X * NUM_FFT_GRID_POINTS_Y * NUM_FFT_GRID_POINTS_Z);

  for (int j = 0; j < fft.local_ny_after_transpose; j++){

    const double q2 = qy(j + fft.local_y_start_after_transpose);

    const int y_index = j*NUM_FFT_GRID_POINTS_X*fft.pad;

    for (int i = 0; i < NUM_FFT_GRID_POINTS_X; i++){

      const double q1 = qx(i);
      const double qq = qysq(j + fft.local_y_start_after_transpose) + qxsq(i) + 1e-16;

      const int x_index = i*fft.pad;

      for (int k = 0; k < fft.pad; k++){

        const int index = x_index + y_index + k;

        const double q3 = qpad(k);
        const double f1 = fx[index];
        const double f2 = fy[index];
        const double f3 = fz[index];
        const double norm = 1.0/(qq + qpadsq(k));
        const double kdotf = (q1*f1+q2*f2+q3*f3)*norm;
        ux[index] = norm*(f1-q1*(kdotf))/(FLUID_VISCOSITY*totalpts);
        uy[index] = norm*(f2-q2*(kdotf))/(FLUID_VISCOSITY*totalpts);
        uz[index] = norm*(f3-q3*(kdotf))/(FLUID_VISCOSITY*totalpts);

      }
    }
  }

  if (fft.local_y_start_after_transpose == 0){

    ux[0] = 0.0;
    ux[1] = 0.0;
    uy[0] = 0.0;
    uy[1] = 0.0;
    uz[0] = 0.0;
    uz[1] = 0.0;

  }

  fft.inverse_transform(ux);
  fft.inverse_transform(uy);
  fft.inverse_transform(uz);

  for (int index = 0; index < fft.total_local_size; index++){

    ux_shear(index) = ux[index];
    uy_shear(index) = uy[index];
    uz_shear(index) = uz[index];

  }

  #endif

}

void fcm_fluid_solver::produce_reference_matrices(collection_of_rigid_bodies& all_bodies){

  #if !SIMULATING_BLOBS

    vec evec(3*all_bodies.total_num_blobs, fill::zeros);

    mat single_body_MFCM(3*all_bodies.blobs_per_body, 3*all_bodies.blobs_per_body);

    all_bodies.fill_all_blob_pos();
    set_up_gaussians(all_bodies.all_blob_pos);

    for (int n=0; n < 3*all_bodies.blobs_per_body; n++){

      evec(std::max<int>(0, n-1)) = 0.0;
      evec(n) = 1.0;

      set_up_force_distribution(evec);
      produce_deterministic_velocity_field();
      calculate_deterministic_blob_velocities();

      single_body_MFCM.col(n) = V.subvec(0, 3*all_bodies.blobs_per_body - 1);

    }

    ref_Minv = inv_sympd(single_body_MFCM);

    mat K(3*all_bodies.blobs_per_body, 6);

    for (int m=0; m < all_bodies.blobs_per_body; m++){

      K(3*m, 0, size(3,3)) = eye(3,3);
      K(3*m, 3, size(3,3)) = rcross(all_bodies.blob_ref_pos.col(m));

    }

    ref_N = inv(K.t() * ref_Minv * K);

  #endif

}
