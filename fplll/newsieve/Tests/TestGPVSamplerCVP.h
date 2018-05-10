#ifndef TEST_GPV_SAMPLER_CVP_H
#define TEST_GPV_SAMPLER_CVP_H


#include "../GPVSamplerCVP.h"

#include "TestSampler.h"

bool test_gpv_sampler_cvp(bool output = true)
{

  using Traits = GaussSieve::DefaultSieveTraits<mpz_class, false, -1>;
  using Sampler = GaussSieve::GPVSamplerCVP<Traits,false,std::mt19937,std::seed_seq>;
  
  return test_general_sampler<Sampler>(output);
}


#endif
