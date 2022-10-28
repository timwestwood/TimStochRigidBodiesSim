// config.hpp

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Include guard
#ifndef MY_CONFIG_HEADER_INCLUDED
#define MY_CONFIG_HEADER_INCLUDED

#define SIMULATION_NAME "example"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Simulation type
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define RIGID_BODY_TYPE 2
// Valid options:
// 0 = blobs
// 1 = tetrahedrons
// 2 = rods
// 3 = czech hedgehogs

#define INTERACTION_POTENTIAL_TYPE 0
// Valid options:
// 0 = Basic soft sphere potential
// 1 = Hydrophobic soft sphere potential
// 2 = Hydrophillic soft sphere potential

#define SIMULATION_ALGORITHM 0
// Valid options:
// 0 = Markov chain Monte Carlo (no dynamics -- samples equilibrium distributions only)
// 1 = Generalised Drifter-Corrector
// 2 = Euler-Maruyama

#define SHEARING_FORCES false // Constant sinusoidal fluid forcing in x direction
#define SLIP_CHANNEL true // Image system in the y direction

#define SKIP_ON_OVERLAP false

#define INITIAL_CONDITIONS_TYPE 2
// Valid options:
// 0 = Use the LAST state in the provided file.
// 1 = Use a RANDOM state in the provided file.
// 2 = Randomly seed the bodies and move according to a drag mobility to avoid overlap.
#define INITIAL_CONDITIONS_FILE_NAME "example_state.dat" // Include file extension!

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Physical parameters
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define NUM_BODIES 3
#define BLOB_RADIUS 1.0
#define FLUID_VISCOSITY 1.0
#define ENERGY_SCALE 1.0 // = k_B * T

#if SHEARING_FORCES

  #define SHEAR_FORCE_SCALE 1.0
  #define SHEAR_WAVENUMBER 1.0

#endif

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Computational parameters
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define GMRES_TOL 1e-3
#define GMRES_MAX_ITER 40

#define NUM_FFT_GRID_POINTS_X 256
#define NUM_FFT_GRID_POINTS_Y 256 // The inhabited y domain will correspond to half this value if SLIP_CHANNEL == true.
#define NUM_FFT_GRID_POINTS_Z 256

#define NUM_FCM_GAUSSIAN_POINTS 20

#define NUM_LINKED_LIST_CELLS_PER_DIM 6

#define FINITE_DIFFERENCE_DELTA_FACTOR 1e-5

#define STEPS_PER_DIFFUSIVE_TIME 100
#define NUM_DIFFUSIVE_TIMES 100

#define DATA_SAVES_PER_DIFFUSIVE_TIME 10

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Derived parameters (these should be left alone)
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#define BLOB_DIFFUSIVE_TIME (PI * FLUID_VISCOSITY * BLOB_RADIUS * BLOB_RADIUS * BLOB_RADIUS / ENERGY_SCALE)
#define DT (BLOB_DIFFUSIVE_TIME / STEPS_PER_DIFFUSIVE_TIME)

#define PI 3.14159265358979323846264338327950288

#define DX (BLOB_RADIUS * 0.3033277330901916) // = BLOB_RADIUS/(1.86 * sqrt(PI))
#define LX (DX * NUM_FFT_GRID_POINTS_X)
#define LY (DX * NUM_FFT_GRID_POINTS_Y)
#define LZ (DX * NUM_FFT_GRID_POINTS_Z)

#define SIMULATING_RODS (RIGID_BODY_TYPE == 2)
#define SIMULATING_TETRAHEDRA (RIGID_BODY_TYPE == 1)
#define SIMULATING_BLOBS (RIGID_BODY_TYPE == 0)
#define SIMULATING_CZECH_HEDGEHOGS (RIGID_BODY_TYPE == 3)

#define SOFT_SPHERE (INTERACTION_POTENTIAL_TYPE == 0)
#define SOFT_SPHERE_HYDROPHOBIC (INTERACTION_POTENTIAL_TYPE == 1)
#define SOFT_SPHERE_HYDROPHILLIC (INTERACTION_POTENTIAL_TYPE == 2)

#define MARKOV_CHAIN_MONTE_CARLO (SIMULATION_ALGORITHM == 0)
#define GEN_DRIFTER_CORRECTOR (SIMULATION_ALGORITHM == 1)
#define EULER_MARUYAMA (SIMULATION_ALGORITHM == 2)

#define IC_READ_IN_LAST (INITIAL_CONDITIONS_TYPE == 0)
#define IC_READ_IN_RANDOM (INITIAL_CONDITIONS_TYPE == 1)
#define IC_SEED_RANDOM (INITIAL_CONDITIONS_TYPE == 2)

#define READ_IN_INITIAL_CONDITIONS (IC_READ_IN_LAST || IC_READ_IN_RANDOM)

#if MARKOV_CHAIN_MONTE_CARLO

  #define DT 1.0

#endif

#endif // MY_CONFIG_HEADER_INCLUDED
