// fft_solver.cpp

#include "fft_solver.hpp"
#include "config.hpp"

fft_solver::~fft_solver(){

  delete[] workspace;

  rfftwnd_mpi_destroy_plan(plan);
  rfftwnd_mpi_destroy_plan(inverse_plan);

}

fft_solver::fft_solver(){}

void fft_solver::initialise(){

  // See FFTW docs (e.g. http://www.fftw.org/fftw2_doc/fftw_4.html) for details

  plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, NUM_FFT_GRID_POINTS_X, NUM_FFT_GRID_POINTS_Y, NUM_FFT_GRID_POINTS_Z, FFTW_REAL_TO_COMPLEX, FFTW_MEASURE);
  inverse_plan = rfftw3d_mpi_create_plan(MPI_COMM_WORLD, NUM_FFT_GRID_POINTS_X, NUM_FFT_GRID_POINTS_Y, NUM_FFT_GRID_POINTS_Z, FFTW_COMPLEX_TO_REAL, FFTW_MEASURE);

  rfftwnd_mpi_local_sizes(plan, &local_nx, &local_x_start, &local_ny_after_transpose, &local_y_start_after_transpose, &total_local_size);

  pad = 2 * (NUM_FFT_GRID_POINTS_Z/2 + 1);

  workspace = new double[total_local_size];

}

void fft_solver::transform(double *in){

  rfftwnd_mpi(plan, 1, in, workspace, FFTW_TRANSPOSED_ORDER);

}

void fft_solver::inverse_transform(double *in){

  rfftwnd_mpi(inverse_plan, 1, in, workspace, FFTW_TRANSPOSED_ORDER);

}
