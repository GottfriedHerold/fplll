#ifndef TEST_QUEUE_H
#define TEST_QUEUE_H

#include <iostream>
#include "fplll.h"
#include "../Typedefs.h"
#include "gmpxx.h"
#include <random>
#include "../LatticeBases.h"
#include "../GaussQueue.h"
#include "../GPVSampler.h"
#include "../GaussQueue_impl.h"

bool test_queue(bool output = true)
{
  int constexpr dim = 20;
  fplll::ZZ_mat<mpz_t> B;
  B.resize(dim, dim);
  srand (1);

  //generates GM lattice
  B.gen_qary_prime(1, 10*dim);

  fplll::lll_reduction(B, fplll::LLL_DEF_DELTA, fplll::LLL_DEF_ETA, fplll::LM_WRAPPER);
  using Traits = GaussSieve::DefaultSieveTraits<mpz_class, false, -1>;
  typename Traits::GlobalStaticDataInitializer init_arg (dim);
  using Sampler = GaussSieve::GPVSampler<Traits,false,std::mt19937,std::seed_seq>;
  GaussSieve::SieveLatticeBasis<Traits,false> sieve_basis(B,init_arg); // convert to SieveLatticeBasis


  // sampler2 is a free-standing sampler for B.
  std::seed_seq sseq {1,2,3,5};
  Sampler sampler2(sseq);
  sampler2.init(nullptr, sieve_basis);

  // gauss_queue contains another sampler.
  GaussSieve::GaussQueue<Traits,false> gauss_queue(nullptr,init_arg, 12345, nullptr);
  // We have to initialize the sampler inside the seed as well.
  gauss_queue.sampler->init(nullptr, sieve_basis);

  if (output) std::cout << "taking 20 elements from the queue (sampling)" << std::endl;
  for(int i=0;i<20;++i)
  {
    auto point = gauss_queue.true_pop();
    if (output) std::cout << point << std::endl;
  }

  if (output) std::cout << "sampling single element" << std::endl;

  auto LP = sampler2.sample(); // LP is non-zero for the above sseq.
  if (output) std::cout << LP << std::endl;
  auto LP2 = LP.make_copy();
  gauss_queue.push(std::move(LP));
  assert(gauss_queue.size()==1);
  if (output) std::cout << "Outputting queue :" << gauss_queue << std::endl;
  LP = gauss_queue.true_pop();
  assert(gauss_queue.size()==0);
  assert(LP==LP2);

  gauss_queue.push(std::move(LP));
  std::stringstream channel;
  if (gauss_queue.serialize_gauss_queue(channel) == false)
  {
    assert(false);
  }
  GaussSieve::GaussQueue<Traits, false> gauss_queue2(channel, nullptr);
  channel.clear();
  if (gauss_queue.serialize_gauss_queue(channel) == false)
  {
    assert(false);
  }
  GaussSieve::GaussQueue<Traits, false> gauss_queue3{nullptr, init_arg, 112233, nullptr};

  gauss_queue3.unserialize_gauss_queue(channel);
  assert(gauss_queue2.size() == 1);
  assert(gauss_queue3.size() == 1);
  auto LP3 = gauss_queue2.true_pop();
  auto LP4 = gauss_queue3.true_pop();
  assert(LP3 == LP4);
  gauss_queue2.sampler->init(nullptr, sieve_basis);
  gauss_queue3.sampler->init(nullptr, sieve_basis);
  LP3 = gauss_queue2.true_pop();
  LP4 = gauss_queue3.true_pop();
  assert(LP3 == LP4);
  return true;
}

#endif
