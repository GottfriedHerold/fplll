#ifndef TEST_BLOCK_ORTHOGONAL_SIMHASH_H
#define TEST_BLOCK_ORTHOGONAL_SIMHASH_H

#include "TestSimHash.h"

#include "../BlockOrthogonalSimHash.h"
#include <sstream>

bool test_block_orthogonal_sim_hash()
{
  using MyCooSelect = GaussSieve::BlockOrthogonalSimHash<64,2,false, long int>;
  // test_general_coo_selection defined in TestSimHash.
  test_general_coo_selection<MyCooSelect>();

//  std::cout << GaussSieve::GlobalBitApproxData<MyCooSelect>::coo_selection;

  std::mt19937 rng; // default-seeded
  GaussSieve::PMatrix pm1{64,rng};
  GaussSieve::PMatrix pm2;
  pm1.print();
  std::stringstream channel;
  channel << pm1;
  channel >> pm2;
  assert(pm1 == pm2);

  GaussSieve::DMatrix dm1{64,rng};
  GaussSieve::DMatrix dm2;
  channel.clear();
  channel << dm1;
  channel >> dm2;
  assert(dm1 == dm2);
  return true;
}

#endif
