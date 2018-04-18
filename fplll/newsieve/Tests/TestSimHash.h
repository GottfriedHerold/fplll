#ifndef TEST_SIMHASHES_H
#define TEST_SIMHASHES_H

#include "../DebugAll.h"
#include "../Typedefs.h"
#include <iostream>
#include <type_traits>
#include "../SieveUtility.h"
#include "../SimHash.h"
#include <limits>
#include "../ExactLatticePoint.h"
#include "gmpxx.h"
//#include <math.h>

bool test_sim_hashes()
{
  static constexpr std::size_t sim_hash_len1 = 24;
  static constexpr std::size_t sim_hash_num1 = 2;
  using SimHashBlock1 = std::bitset<sim_hash_len1>; // intentionally not a power of 2.
  using SimHash1      = std::array<SimHashBlock1,sim_hash_num1>;
  SimHash1 sim_hash1;
  for(auto i =0; i<24; ++i)
  {
    sim_hash1[0][i] = 0;
    sim_hash1[1][i] = 1;
  }
  sim_hash1[0][9]  = 1;
  sim_hash1[0][10] = 1;
  sim_hash1[1][1]  = 0;
  assert(sim_hash1[0].count() == 2);
  assert(sim_hash1[1].count() == 23);
  SimHash1 sim_hash2 = GaussSieve::flip_all_bits(sim_hash1);
  assert(sim_hash2[0].count() == 22);
  assert(sim_hash2[1].count() == 1);
  using GaussSieve::operator<<;
  std::cout << "SimHash1:" << sim_hash1 << std::endl;

  std::cout << "Negated :" << sim_hash2 << std::endl;
  return true;
}

// test for general property of any coo_selection
namespace TestHelper
{
  // these have to be at NS level rather than inside the function where we use them.
  template<class PurportedCooSelection> using GetSimHashBlock  = typename PurportedCooSelection::SimHashBlock;
  template<class PurportedCooSelection> using GetSimHashes     = typename PurportedCooSelection::SimHashes;
  template<class PurportedCooSelection> using GetDimensionType = typename PurportedCooSelection::DimensionType;
} // end namespace TestHelper

template<class CoordinateSelection>
bool test_general_coo_selection()
{
  // static checks:

  static_assert(GaussSieve::IsACoordinateSelection<CoordinateSelection>::value,"");
  // checks for the existence of typedefs SimHashBlock, SimHashes in CooSelection.
  static_assert(GaussSieve::mystd::is_detected<TestHelper::GetSimHashBlock,  CoordinateSelection>::value,"");
  static_assert(GaussSieve::mystd::is_detected<TestHelper::GetSimHashes,     CoordinateSelection>::value,"");
  static_assert(GaussSieve::mystd::is_detected<TestHelper::GetDimensionType, CoordinateSelection>::value,"");

  using SimHashBlock = typename CoordinateSelection::SimHashBlock;
  using SimHashes    = typename CoordinateSelection::SimHashes;
  using DimensionType= typename CoordinateSelection::DimensionType;
  using SimHashDeref = GaussSieve::mystd::decay_t<decltype (std::declval<SimHashes>()[0])>;
  static_assert(std::is_same<SimHashBlock, SimHashDeref>::value, "");
  static_assert(std::is_default_constructible<CoordinateSelection>::value,"");
  static_assert(std::is_copy_assignable<CoordinateSelection>::value,"");
  static_assert(std::is_default_constructible<typename CoordinateSelection::DimensionType>::value,"");
  static_assert(std::is_copy_assignable<typename CoordinateSelection::DimensionType>::value,"");

  // create 2 lattice point latp, latp2:

  int constexpr dim = 67;
  using LP = GaussSieve::ExactLatticePoint<long, -1>;
  using GaussSieve::MaybeFixed;
  std::array<long,dim> arr;
  std::array<long,dim> arr2;
  for(int i=0;i<dim;++i)
  {
    arr[i] = std::pow(-1, i+1) * i;
    arr2[i] = std::pow(-1, i) * (i-1)+13;
  }

  GaussSieve::StaticInitializer<LP> init1 (MaybeFixed<-1>{dim});
  LP latp = GaussSieve::make_from_any_vector<LP>(arr,MaybeFixed<-1>{dim});
  LP latp2 = GaussSieve::make_from_any_vector<LP>(arr2,MaybeFixed<-1>{dim});
  std::cout << "latp = " << latp << std::endl;
  std::cout << "latp2 = " << latp2 << std::endl;

  // initialize the hash function with some arbitrary random seed.
  GaussSieve::StaticInitializer<GaussSieve::GlobalBitApproxData<CoordinateSelection>> init_coo_selection(DimensionType{dim},5235);

  SimHashes sh1 = GaussSieve::GlobalBitApproxData<CoordinateSelection>::coo_selection.compute_all_bitapproximations(latp);
  SimHashes sh2 = GaussSieve::GlobalBitApproxData<CoordinateSelection>::coo_selection.compute_all_bitapproximations(latp2);
  using GaussSieve::operator<<;
  std::cout << "SimHash for latp = " << sh1 << std::endl;
  std::cout << "SimHash for latp2= " << sh2 << std::endl;

  using LPWithBitApprox = GaussSieve::AddBitApproximationToLP<LP, CoordinateSelection>;
  LPWithBitApprox latp_with_bitapprox2{std::move(latp2)};
  std::cout << "Point2 with Bitapprox " << latp_with_bitapprox2 << std::endl;

  return true;
}

#endif
