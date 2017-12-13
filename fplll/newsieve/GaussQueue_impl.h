#ifndef GAUSS_QUEUE_CPP
#define GAUSS_QUEUE_CPP

#include "GaussQueue.h"
#include "GPVSampler.h"
#include "UniformSampler.h"
#include "GPVSamplerExtended.h"
#include "Sampler.h"
#include "DebugAll.h"

namespace GaussSieve
{

// constructor, single-threaded version
template<class SieveTraits>
GaussQueue<SieveTraits,false>::GaussQueue(Sieve<SieveTraits,false> * const caller_sieve,
                                          GlobalStaticDataInitializer const &static_data,
                                          int seed_sampler,
                                          Sampler<SieveTraits,false,std::mt19937_64,std::seed_seq> *user_sampler)
    : init_data_type(static_data),
      init_ret_type(static_data),
      main_queue(),
      gauss_sieve(caller_sieve),
      sampler(user_sampler),
      sampler_owned(user_sampler == nullptr)
{
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
  assert(caller_sieve==nullptr);
#else
  assert(caller_sieve!=nullptr);
#endif
  if (sampler == nullptr)
  {
    // Use seed_sampler to init the rng. The other numbers spell "SAMPLER" in ASCII, they are here
    //  because we use the seed elsewhere as well.
    std::seed_seq seed{83, 65, 77, 80, 76, 69, 82, seed_sampler};
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing our own sampler. Using GPVSampler")
    sampler = new GPVSampler<SieveTraits, false, std::mt19937_64, std::seed_seq>(seed);
  }
  assert(sampler != nullptr);
}

template<class SieveTraits>
auto GaussQueue<SieveTraits,false>::true_pop() -> RetType
{
  if(main_queue.empty()) // Queue is empty, sample a new element.
  {
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
    assert(gauss_sieve==nullptr);
#else
    assert(gauss_sieve!=nullptr);
    gauss_sieve->statistics.increment_number_of_points_sampled();
    gauss_sieve->statistics.increment_number_of_points_constructed();
#endif
    assert(sampler!=nullptr);
    RetType ret{sampler->sample()};
    return ret;
    //return static_cast<typename SieveTraits::GaussList_StoredPoint>(sampler->sample());
  }
  else  // Queue is not empty, just return "next" stored element.
        // (for std::priority_queue, this is called front(), for normal std::queue(), it's top().
  {
#ifndef USE_REGULAR_QUEUE
    RetType ret{std::move(*(main_queue.top())};  // move from the queue.
    // Note: The top of the queue still holds a valid pointer to a lattice point
    // the std::move above just put that lattice point into an unspecified and unusable state.
    // we still need to free its memory.
    delete main_queue.top();
#else   // analogous
    RetType ret{std::move(*(main_queue.front()))};
    delete main_queue.front();
#endif // USE_REGULAR_QUEUE
    main_queue.pop();  //This removes the pointer itself from the queue.
    return ret;
  }
}

template<class SieveTraits>
void GaussQueue<SieveTraits,false>::push(DataType && val)
{
  DataType * new_lp = new DataType (std::move(val) );
  main_queue.push(new_lp);
}



template<class SieveTraits> GaussQueue<SieveTraits,false>::~GaussQueue()
{
  while(! main_queue.empty() )
  {
#ifndef USE_REGULAR_QUEUE
  delete main_queue.top();
#else
  delete main_queue.front();
#endif // USE_REGULAR_QUEUE
  main_queue.pop();
  }
  if (sampler_owned)
  {
    delete sampler;
  }
}

}  // end namespace GaussSieve

#endif  // include guard
