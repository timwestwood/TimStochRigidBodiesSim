// linked_list.cpp

#include <cmath>
#include <armadillo>
#include "linked_list.hpp"
#include "config.hpp"
using namespace arma;

linked_list::~linked_list(){}

linked_list::linked_list(){}

int linked_list::get_cell_id(const int i, const int j, const int k) const {

  const int imod = int(fmod(i+NUM_LINKED_LIST_CELLS_PER_DIM,NUM_LINKED_LIST_CELLS_PER_DIM));
  const int jmod = int(fmod(j+NUM_LINKED_LIST_CELLS_PER_DIM,NUM_LINKED_LIST_CELLS_PER_DIM));
  const int kmod = int(fmod(k+NUM_LINKED_LIST_CELLS_PER_DIM,NUM_LINKED_LIST_CELLS_PER_DIM));

  return imod + jmod*NUM_LINKED_LIST_CELLS_PER_DIM + kmod*NUM_LINKED_LIST_CELLS_PER_DIM*NUM_LINKED_LIST_CELLS_PER_DIM;

}

void linked_list::initialise(const int total_num_blobs, const int myrank, const int totalnodes){

  const int total_num_cells = NUM_LINKED_LIST_CELLS_PER_DIM*NUM_LINKED_LIST_CELLS_PER_DIM*NUM_LINKED_LIST_CELLS_PER_DIM;

  my_start = myrank * int(double(total_num_cells)/double(totalnodes));
  my_end = (1 + myrank) * int(double(total_num_cells)/double(totalnodes));

  if (myrank == totalnodes-1){

    my_end = total_num_cells;

  }

  head = std::vector<int>(total_num_cells);

  list = std::vector<int>(total_num_blobs);

  map = std::vector<int>(13*total_num_cells);

  for (int i=0; i<NUM_LINKED_LIST_CELLS_PER_DIM; i++){

    for (int j=0; j<NUM_LINKED_LIST_CELLS_PER_DIM; j++){

      for (int k=0; k<NUM_LINKED_LIST_CELLS_PER_DIM; k++){

        int my_cell = get_cell_id(i,j,k);

        int map_pos = 13*my_cell;

        map[map_pos] = get_cell_id(i+1,j,k);
        map[map_pos+1] = get_cell_id(i+1,j+1,k);
        map[map_pos+2] = get_cell_id(i,j+1,k);
        map[map_pos+3] = get_cell_id(i-1,j+1,k);
        map[map_pos+4] = get_cell_id(i+1,j,k-1);
        map[map_pos+5] = get_cell_id(i+1,j+1,k-1);
        map[map_pos+6] = get_cell_id(i,j+1,k-1);
        map[map_pos+7] = get_cell_id(i-1,j+1,k-1);
        map[map_pos+8] = get_cell_id(i+1,j,k+1);
        map[map_pos+9] = get_cell_id(i+1,j+1,k+1);
        map[map_pos+10] = get_cell_id(i,j+1,k+1);
        map[map_pos+11] = get_cell_id(i-1,j+1,k+1);
        map[map_pos+12] = get_cell_id(i,j,k+1);

      }

    }

  }

}

int linked_list::find_containing_cell(const mat& Y) const {

  return floor(NUM_LINKED_LIST_CELLS_PER_DIM*(Y(0)/LX)) + floor(NUM_LINKED_LIST_CELLS_PER_DIM*(Y(1)/LY))*NUM_LINKED_LIST_CELLS_PER_DIM + floor(NUM_LINKED_LIST_CELLS_PER_DIM*(Y(2)/LZ))*NUM_LINKED_LIST_CELLS_PER_DIM*NUM_LINKED_LIST_CELLS_PER_DIM;

}

void linked_list::populate(const mat& Y){

  const int total_num_cells = NUM_LINKED_LIST_CELLS_PER_DIM*NUM_LINKED_LIST_CELLS_PER_DIM*NUM_LINKED_LIST_CELLS_PER_DIM;

  for (int i=0; i<total_num_cells; i++){

    head[i] = -1;

  }

  const int total_num_blobs = list.size();

  for (int i=0; i<total_num_blobs; i++){

    int my_cell = find_containing_cell(Y.row(i));

    list[i] = head[my_cell];

    head[my_cell] = i;

  }

}

bool linked_list::check_for_overlap(const mat& Y) const {

  const int blobs_per_body = list.size()/NUM_BODIES;

  for (int i=my_start; i<my_end; i++){

    int blob1 = head[i];

    while (blob1 != -1){

      const double x1 = Y(blob1,0);
      const double y1 = Y(blob1,1);
      const double z1 = Y(blob1,2);

      const int body1 = blob1/blobs_per_body;

      // same cell
      int blob2 = list[blob1];

      while (blob2 != -1){

        const int body2 = blob2/blobs_per_body;

        if (body1 != body2){

          const double x2 = Y(blob2,0);
          const double y2 = Y(blob2,1);
          const double z2 = Y(blob2,2);

          const double dist2 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);

          if (dist2 < 4.0*BLOB_RADIUS*BLOB_RADIUS){

            return true;

          }

        }

        blob2 = list[blob2];

      }

      // neighbour cells
      for (int j=0; j<13; j++){

        blob2 = head[map[13*i + j]];

        while (blob2 != -1){

          const int body2 = blob2/blobs_per_body;

          if (body1 != body2){

            const double x2 = Y(blob2,0);
            const double y2 = Y(blob2,1);
            const double z2 = Y(blob2,2);

            const double xdiff = x1 - x2 - LX*trunc(2.0*(x1 - x2)/LX);
            const double ydiff = y1 - y2 - LY*trunc(2.0*(y1 - y2)/LY);
            const double zdiff = z1 - z2 - LZ*trunc(2.0*(z1 - z2)/LZ);

            const double dist2 = xdiff*xdiff + ydiff*ydiff + zdiff*zdiff;

            if (dist2 < 4.0*BLOB_RADIUS*BLOB_RADIUS){

              return true;

            }

          }

          blob2 = list[blob2];

        }

      }

      blob1 = list[blob1];

    }

  }

  return false;

}

void linked_list::apply_barrier_forces(vec& blob_forces, const mat& Y) const {

  const int blobs_per_body = list.size()/NUM_BODIES;

  for (int i=my_start; i<my_end; i++){

    int blob1 = head[i];

    while (blob1 != -1){

      const double x1 = Y(blob1,0);
      const double y1 = Y(blob1,1);
      const double z1 = Y(blob1,2);

      const int body1 = blob1/blobs_per_body;

      // same cell
      int blob2 = list[blob1];

      while (blob2 != -1){

        const int body2 = blob2/blobs_per_body;

        if (body1 != body2){ // Don't repel blobs in the same rigid body

          const double x2 = Y(blob2,0);
          const double y2 = Y(blob2,1);
          const double z2 = Y(blob2,2);

          const double r = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

          const double magnitude = interaction_magnitude(r);

          const vec3 direction = {(x1-x2)/r, (y1-y2)/r, (z1-z2)/r};

          blob_forces.subvec(3*blob1, 3*blob1 + 2) += magnitude*direction;

          blob_forces.subvec(3*blob2, 3*blob2 + 2) -= magnitude*direction;

        }

        blob2 = list[blob2];

      }

      // neighbour cells
      for (int j=0; j<13; j++){

        blob2 = head[map[13*i + j]];

        while (blob2 != -1){

          const int body2 = blob2/blobs_per_body;

          if (body1 != body2){ // Don't repel blobs in the same rigid body

            const double x2 = Y(blob2,0);
            const double y2 = Y(blob2,1);
            const double z2 = Y(blob2,2);

            // nearest image
            const double xdiff = x1 - x2 - LX*trunc(2.0*(x1 - x2)/LX);
            const double ydiff = y1 - y2 - LY*trunc(2.0*(y1 - y2)/LY);
            const double zdiff = z1 - z2 - LZ*trunc(2.0*(z1 - z2)/LZ);

            const double r = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

            const double magnitude = interaction_magnitude(r);

            const vec3 direction = {xdiff/r, ydiff/r, zdiff/r};

            blob_forces.subvec(3*blob1, 3*blob1 + 2) += magnitude*direction;

            blob_forces.subvec(3*blob2, 3*blob2 + 2) -= magnitude*direction;

          }

          blob2 = list[blob2];

        }

      }

      blob1 = list[blob1];

    }

  }

}

void linked_list::apply_barrier_potential(double& potential, const arma::mat& Y) const {

  const int blobs_per_body = list.size()/NUM_BODIES;

  for (int i=my_start; i<my_end; i++){

    int blob1 = head[i];

    while (blob1 != -1){

      const double x1 = Y(blob1,0);
      const double y1 = Y(blob1,1);
      const double z1 = Y(blob1,2);

      const int body1 = blob1/blobs_per_body;

      // same cell
      int blob2 = list[blob1];

      while (blob2 != -1){

        const int body2 = blob2/blobs_per_body;

        if (body1 != body2){ // Don't repel blobs in the same rigid body

          const double x2 = Y(blob2,0);
          const double y2 = Y(blob2,1);
          const double z2 = Y(blob2,2);

          const double r = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2));

          potential += interaction_magnitude(r);

        }

        blob2 = list[blob2];

      }

      // neighbour cells
      for (int j=0; j<13; j++){

        blob2 = head[map[13*i + j]];

        while (blob2 != -1){

          const int body2 = blob2/blobs_per_body;

          if (body1 != body2){ // Don't repel blobs in the same rigid body

            const double x2 = Y(blob2,0);
            const double y2 = Y(blob2,1);
            const double z2 = Y(blob2,2);

            // nearest image
            const double xdiff = x1 - x2 - LX*trunc(2.0*(x1 - x2)/LX);
            const double ydiff = y1 - y2 - LY*trunc(2.0*(y1 - y2)/LY);
            const double zdiff = z1 - z2 - LZ*trunc(2.0*(z1 - z2)/LZ);

            const double r = sqrt(xdiff*xdiff + ydiff*ydiff + zdiff*zdiff);

            potential += interaction_magnitude(r);

          }

          blob2 = list[blob2];

        }

      }

      blob1 = list[blob1];

    }

  }

}

double linked_list::interaction_magnitude(const double r) const {

  double magnitude = 0.0;
  // The output will be be value of the potential in MCMC simulations and the magnitude of the force otherwise.

  #if SOFT_SPHERE

    // This potential has the form U(r) = U0 + U0*(2a - r)/b when particles overlap and
    // U(r) = U0*exp((2a - r)/b) otherwise.
    const double U0 = 20.0*ENERGY_SCALE;
    const double b = 0.5*BLOB_RADIUS;

    #if MARKOV_CHAIN_MONTE_CARLO

      magnitude = (r < 2.0*BLOB_RADIUS) ? U0 + U0*(2.0*BLOB_RADIUS - r)/b : U0*exp((2.0*BLOB_RADIUS - r)/b);

    #else

      magnitude = (r < 2.0*BLOB_RADIUS) ? U0/b : U0*exp((2.0*BLOB_RADIUS - r)/b)/b;

    #endif

  #elif SOFT_SPHERE_HYDROPHOBIC

    // This potential has the form U(r) = U0*exp((2a - r)/b) - A*exp(-((r-3a)/lambda)^2)
    // for non-overlapping particles, and is extended as U(2a) - U'(2a)*(2a - r) for r < 2a so
    // that the force magntiude remains constant when particles overlap.
    const double U0 = 20.0*ENERGY_SCALE;
    const double b = 0.5*BLOB_RADIUS;
    const double A = 5.0*ENERGY_SCALE;
    const double lambda = 0.5*BLOB_RADIUS;

    #if MARKOV_CHAIN_MONTE_CARLO

      magnitude = (r < 2.0*BLOB_RADIUS) ? U0 + U0*(2.0*BLOB_RADIUS - r)/b : U0*exp((2.0*BLOB_RADIUS - r)/b);

      if (r < 2.0*BLOB_RADIUS){

        double temp = BLOB_RADIUS/lambda;
        temp = A*exp(-temp*temp);
        magnitude -= temp; // Correct the U(2a) part.
        magnitude += 2.0*BLOB_RADIUS*(2.0*BLOB_RADIUS - r)*temp/(lambda*lambda); // Correct the -U'(2a)*(2a - r) part.

      } else {

        const double temp = (r - 3.0*BLOB_RADIUS)/lambda;
        magnitude -= A*exp(-temp*temp);

      }

    #else

      magnitude = (r < 2.0*BLOB_RADIUS) ? U0/b : U0*exp((2.0*BLOB_RADIUS - r)/b)/b;

      if (r < 2.0*BLOB_RADIUS){

        double temp = BLOB_RADIUS/lambda;
        temp = A*exp(-temp*temp);
        magnitude += 2.0*BLOB_RADIUS*temp/(lambda*lambda);

      } else {

        double temp = (r - 3.0*BLOB_RADIUS)/lambda;
        magnitude += -2.0*A*temp*exp(-temp*temp)/lambda;

      }

    #endif

  #elif SOFT_SPHERE_HYDROPHILLIC

    // This potential has the form U(r) = U0*exp((2a - r)/b) - A*exp(-((r-3a)/lambda)^2) + A*exp(-((r - 7a/2)/lambda)^2)
    // for non-overlapping particles, and is extended as U(2a) - U'(2a)*(2a - r) for r < 2a so
    // that the force magntiude remains constant when particles overlap.
    const double U0 = 20.0*ENERGY_SCALE;
    const double b = 0.5*BLOB_RADIUS;
    const double A = 5.0*ENERGY_SCALE;
    const double lambda = 0.5*BLOB_RADIUS;

    #if MARKOV_CHAIN_MONTE_CARLO

      magnitude = (r < 2.0*BLOB_RADIUS) ? U0 + U0*(2.0*BLOB_RADIUS - r)/b : U0*exp((2.0*BLOB_RADIUS - r)/b);

      if (r < 2.0*BLOB_RADIUS){

        double temp = BLOB_RADIUS/lambda;
        temp = A*exp(-temp*temp);
        magnitude -= temp; // Correct the U(2a) part.
        magnitude += 2.0*BLOB_RADIUS*(2.0*BLOB_RADIUS - r)*temp/(lambda*lambda); // Correct the -U'(2a)*(2a - r) part.

        temp = 1.5*BLOB_RADIUS/lambda;
        temp = A*exp(-temp*temp);
        magnitude += temp; // Correct the U(2a) part.
        magnitude -= 3.0*BLOB_RADIUS*(2.0*BLOB_RADIUS - r)*temp/(lambda*lambda); // Correct the -U'(2a)*(2a - r) part.

      } else {

        double temp = (r - 3.0*BLOB_RADIUS)/lambda;
        magnitude -= A*exp(-temp*temp);

        temp = (r - 3.5*BLOB_RADIUS)/lambda;
        magnitude += A*exp(-temp*temp);

      }

    #else

      magnitude = (r < 2.0*BLOB_RADIUS) ? U0/b : U0*exp((2.0*BLOB_RADIUS - r)/b)/b;

      if (r < 2.0*BLOB_RADIUS){

        double temp = BLOB_RADIUS/lambda;
        temp = A*exp(-temp*temp);
        magnitude += 2.0*BLOB_RADIUS*temp/(lambda*lambda);

        temp = 1.5*BLOB_RADIUS/lambda;
        temp = A*exp(-temp*temp);
        magnitude -= 3.0*BLOB_RADIUS*temp/(lambda*lambda);

      } else {

        double temp = (r - 3.0*BLOB_RADIUS)/lambda;
        magnitude += -2.0*A*temp*exp(-temp*temp)/lambda;

        temp = (r - 3.5*BLOB_RADIUS)/lambda;
        magnitude += 2.0*A*temp*exp(-temp*temp)/lambda;

      }

    #endif

  #endif

  return magnitude;

}
