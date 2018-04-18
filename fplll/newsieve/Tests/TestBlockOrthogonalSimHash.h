#ifndef TEST_BLOCK_ORTHOGONAL_SIMHASH_H
#define TEST_BLOCK_ORTHOGONAL_SIMHASH_H

#include "TestSimHash.h"

#include "../BlockOrthogonalSimHash.h"

bool test_block_orthogonal_sim_hash()
{
  using MyCooSelect = GaussSieve::BlockOrthogonalSimHash<64,2,false, long int>;
  // test_general_coo_selection defined in TestSimHash.
  return test_general_coo_selection<MyCooSelect>();
}

#endif
