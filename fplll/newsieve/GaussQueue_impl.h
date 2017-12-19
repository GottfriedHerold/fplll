#ifndef GAUSS_QUEUE_CPP
#define GAUSS_QUEUE_CPP

#include "DebugAll.h"
#include "GPVSampler.h"
#include "GPVSamplerExtended.h"
#include "GaussQueue.h"
#include "Sampler.h"
#include "UniformSampler.h"
#include "GaussVectorBitApprox.h"

namespace GaussSieve
{

// constructor, single-threaded version
// clang-format off
template <class SieveTraits>
GaussQueue<SieveTraits, false>::GaussQueue(Sieve<SieveTraits, false> *const caller_sieve,
                                          GaussVectorWithBitApprox<SieveTraits, false> * const caller_list,
                                          GlobalStaticDataInitializer const &static_data,
                                          int seed_sampler,
                                          Sampler<SieveTraits, false, std::mt19937_64, std::seed_seq> *user_sampler)
    : // init_data_type(static_data),
      // init_ret_type(static_data),
      main_queue(),
      gauss_sieve(caller_sieve),
      main_list(caller_list),
      sampler(user_sampler),
      sampler_owned(user_sampler == nullptr)
// clang-format on
{
  assert(caller_list != nullptr);
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
  assert(caller_sieve == nullptr);
#else
  assert(gauss_sieve != nullptr);
  assert(& (gauss_sieve->main_vector) == main_list );
#endif
  if (sampler == nullptr)
  {
    // Use seed_sampler to init the rng. The other numbers spell "SAMPLER" in ASCII, they are here
    //  because we use the very same seed elsewhere as well (Think of sees_seq as a Hash)
    std::seed_seq seed{83, 65, 77, 80, 76, 69, 82, seed_sampler};
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing our own sampler. Using GPVSamplerExtended")
    sampler = new GPVSamplerExtended<SieveTraits, false, std::mt19937_64, std::seed_seq>(seed, caller_sieve->get_ambient_dimension() / 3 );
  }
  assert(sampler != nullptr);
}

template <class SieveTraits> auto GaussQueue<SieveTraits, false>::true_pop() -> StoredIterator
{
  if (main_queue.empty())  // Queue is empty, sample a new element.
  {
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
    assert(gauss_sieve == nullptr);
#else
    assert(gauss_sieve != nullptr);
    gauss_sieve->statistics.increment_number_of_points_sampled();
    gauss_sieve->statistics.increment_number_of_points_constructed();
#endif
    assert(sampler != nullptr);

    auto new_point = sampler->sample();
    return main_list.emplace_back(std::move(new_point));


//    RetType ret{sampler->sample()};
//    return ret;
    // return static_cast<typename SieveTraits::GaussList_StoredPoint>(sampler->sample());
  }
  else  // Queue is not empty, just return "next" stored element.
        // (for std::priority_queue, this is called front(), for normal std::queue(), it's top().
  {
// clang-format off
#ifndef USE_REGULAR_QUEUE
#error update
    RetType ret{std::move(*(main_queue.top()))};  // move from the queue.
    // Note: The top of the queue still holds a valid pointer to a lattice point
    // the std::move above just put that lattice point into an unspecified and unusable state.
    // we still need to free its memory.
    delete main_queue.top();
#else
    StoredIterator ret{std::move(main_queue.front())};
#endif  // USE_REGULAR_QUEUE
// clang-format off
    main_queue.pop();  // This just removes the pointer.
    return ret;
  }
}

/*
// clang-format off
template<class SieveTraits>
void GaussQueue<SieveTraits, false>::push(DataType &&val)
{
#ifndef USE_REGULAR_QUEUE
  DataType *new_lp_ptr = new DataType(std::move(val));
  main_queue.push(new_lp_ptr);
#else
  main_queue.push(std::move(val));
#endif
}
*/

template<class SieveTraits>
inline void GaussQueue<SieveTraits, false>::push(StoredIterator const &val)
{
  main_queue.push(val);
}

template <class SieveTraits>
GaussQueue<SieveTraits, false>::~GaussQueue()
{
// free memory if the queue stores pointers
#ifndef USE_REGULAR_QUEUE
#error Update that to vector as main storage
  while (!main_queue.empty())
  {
    delete main_queue.top();
    main_queue.pop();
  }
#endif
  if (sampler_owned)
  {
    delete sampler;
  }
}
// clang-format on

}  // end namespace GaussSieve

#endif  // include guard
