#ifndef TEST_SAMPLER_H
#define TEST_SAMPLER_H

#include "../DefaultIncludes.h"
#include "../Compat.h"
#include <random>
#include "../LatticeBases.h"
#include <sstream>
#include <iostream>
#include "fplll.h"
#include "../Typedefs.h"
#include "gmpxx.h"

namespace GaussSieveTestHelper
{
template<class T> using GetParent = typename T::Parent;
}


template<class SamplerClass>
bool test_general_sampler(bool output = true)
{

  static_assert(GaussSieve::mystd::is_detected<GaussSieveTestHelper::GetParent, SamplerClass>::value, "Derived sampler must declare a Parent typedef");

  using EngineUsed = typename SamplerClass::Parent::GetEngine;
  using SseqUsed = typename SamplerClass::Parent::GetSseq;
  static_assert(std::is_same<EngineUsed, std::mt19937>::value, "Testing only for default Engine");
  static_assert(std::is_same<SseqUsed, std::seed_seq>::value, "Testing only for default Sseq");
  using Traits = typename SamplerClass::Parent::GetTraits;

  int constexpr dim = 20;
  fplll::ZZ_mat<mpz_t> B;
  B.resize(dim, dim);
  srand (1);

  //generates GM lattice
  B.gen_qary_prime(1, 10*dim);

  fplll::lll_reduction(B, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::LM_WRAPPER);

  if (output)
  {
    std::cout << "Generating random basis: Result is" << std::endl;
    std::cout << B << std::endl;
  }

  typename Traits::GlobalStaticDataInitializer init_arg(dim);
  GaussSieve::SieveLatticeBasis<Traits,false> sieve_basis(B, init_arg); // convert to SieveLatticeBasis

  std::seed_seq sseq {1,2,3};

  SamplerClass sampler(sseq);
  sampler.init(nullptr, sieve_basis);
  for(int i=0; i<100 ; ++i)
  {
    auto res = sampler.sample();
    if(output)
    {
      std::cout << res << std::endl;
    }
  }

  std::stringstream channel;
  channel << sampler;
  SamplerClass sampler2(channel);
  sampler2.init(nullptr, sieve_basis);
  auto res = sampler2.sample();
  if(output)
  {
    std::cout << res << std::endl;
  }
  channel.clear();
  channel << sampler;
  SamplerClass sampler3(sseq);
  channel >> sampler3;
  sampler3.init(nullptr, sieve_basis);
  auto res2 = sampler3.sample();
  if(output)
  {
    std::cout << res2 << std::endl;
  }
  assert(res == res2);
  return true;
}

#endif
