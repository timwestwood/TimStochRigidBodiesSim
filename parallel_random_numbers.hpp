// parallel_random_numbers.hpp

// =============================================================================
// Include guard
#ifndef MY_PARALLEL_RANDOM_NUMBERS_HEADER_INCLUDED
#define MY_PARALLEL_RANDOM_NUMBERS_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
extern "C" {
  #define SIMPLE_SPRNG /* simple interface */
  #define USE_MPI /* SPRNG makes MPI calls */
  #include "sprng.h" /* SPRNG header file */
}

double gaussian_random();

#endif // MY_PARALLEL_RANDOM_NUMBERS_HEADER_INCLUDED
