// fft_solver.hpp

// =============================================================================
// Include guard
#ifndef MY_FFT_SOLVER_HEADER_INCLUDED
#define MY_FFT_SOLVER_HEADER_INCLUDED

#include <rfftw_mpi.h>

class fft_solver {

public:

  rfftwnd_mpi_plan plan;
  rfftwnd_mpi_plan inverse_plan;
  double *workspace;
  int local_nx;
  int local_x_start;
  int local_ny_after_transpose;
  int local_y_start_after_transpose;
  int total_local_size;
  int pad;

  ~fft_solver();
  fft_solver();

  void initialise();
  void transform(double *in);
  void inverse_transform(double *in);

}; // End of class

#endif // MY_FFT_SOLVER_HEADER_INCLUDED
