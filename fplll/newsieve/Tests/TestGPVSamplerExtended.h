#ifndef TEST_GPV_SAMPLER_EXTENDED_H
#define TEST_GPV_SAMPLER_EXTENDED_H


#include "../GPVSamplerExtended.h"

#include "TestSampler.h"

bool test_gpv_sampler_extended(bool output = true)
{

  using Traits = GaussSieve::DefaultSieveTraits<mpz_class, false, -1>;
  using Sampler = GaussSieve::GPVSamplerExtended<Traits,false,std::mt19937,std::seed_seq>;
  // Note: This uses start_babai == 0, (i.e. does not actually use start_babai)
  return test_general_sampler<Sampler>(output);
}


#endif
