// linked_list.hpp

// =============================================================================
// Include guard
#ifndef MY_LINKED_LIST_HEADER_INCLUDED
#define MY_LINKED_LIST_HEADER_INCLUDED

// =============================================================================
// Forward declared dependencies

// =============================================================================
// Included dependencies
#include <stdlib.h>

class linked_list {

public:

  std::vector<int> head;
  std::vector<int> list;
  std::vector<int> map;
  int my_start;
  int my_end;

  ~linked_list();
  linked_list();

  int get_cell_id(const int i, const int j, const int k) const;
  void initialise(const int total_num_blobs, const int myrank, const int totalnodes);
  int find_containing_cell(const arma::mat& Y) const;
  void populate(const arma::mat& Y);
  bool check_for_overlap(const arma::mat& Y) const;
  void apply_barrier_forces(arma::vec& blob_forces, const arma::mat& Y) const;
  void apply_barrier_potential(double& potential, const arma::mat& Y) const;
  double interaction_magnitude(const double r) const;

}; // End of class

#endif // MY_LINKED_LIST_HEADER_INCLUDED
