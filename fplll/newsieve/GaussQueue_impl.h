#ifndef GAUSS_QUEUE_CPP
#define GAUSS_QUEUE_CPP

#include "DebugAll.h"
#include "GPVSampler.h"
#include "GPVSamplerExtended.h"
#include "GaussQueue.h"
#include "Sampler.h"
#include "UniformSampler.h"

namespace GaussSieve
{

// constructor, single-threaded version
// clang-format off
template <class SieveTraits>
GaussQueue<SieveTraits, false>::GaussQueue(Sieve<SieveTraits, false> *const caller_sieve,
                                          GlobalStaticDataInitializer const &static_data,
                                          int new_seed_sampler,
                                          Sampler<SieveTraits, false, std::mt19937_64, std::seed_seq> *user_sampler)
    : init_data_type(static_data),
      init_ret_type(static_data),
      main_queue(),
      gauss_sieve(caller_sieve),
      seed_sampler(new_seed_sampler),
      sampler(user_sampler),
      sampler_auto_created(false)
// clang-format on
{
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
  assert(caller_sieve == nullptr);
#else
  assert(caller_sieve != nullptr);
#endif
  if (sampler == nullptr)
  {
    auto_gen_sampler(); // create our own sampler and set sampler_auto_created
  }
  assert(sampler != nullptr);
}

template <class SieveTraits> auto GaussQueue<SieveTraits, false>::true_pop() -> RetType
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
    RetType ret{sampler->sample()};
    return ret;
    // return static_cast<typename SieveTraits::GaussList_StoredPoint>(sampler->sample());
  }
  else  // Queue is not empty, just return "next" stored element.
        // (for std::priority_queue, this is called front(), for normal std::queue(), it's top().
  {
// clang-format off
#ifndef USE_REGULAR_QUEUE
    RetType ret{std::move(*(main_queue.top()))};  // move from the queue.
    // Note: The top of the queue still holds a valid pointer to a lattice point
    // the std::move above just put that lattice point into an unspecified and unusable state.
    // we still need to free its memory.
    delete main_queue.top();
#else
    RetType ret{std::move(main_queue.front())};
#endif  // USE_REGULAR_QUEUE
// clang-format off
    main_queue.pop();  // This just removes the pointer.
    return ret;
  }
}

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

template <class SieveTraits>
GaussQueue<SieveTraits, false>::~GaussQueue()
{
// free memory if the queue stores pointers
#ifndef USE_REGULAR_QUEUE
  while (!main_queue.empty())
  {
    delete main_queue.top();
    main_queue.pop();
  }
#endif
  if (sampler_auto_created)
  {
    delete sampler;
  }
}

template<class SieveTraits>
auto GaussQueue<SieveTraits, false>::set_sampler(SamplerPtr const &new_sampler_ptr) -> SamplerPtr
{
  if(sampler_auto_created)
  {
    delete sampler;
    sampler = nullptr;
    sampler_auto_created = false;
  }
  SamplerPtr old_val = sampler;
  sampler = new_sampler_ptr;
  if(sampler == nullptr)
  {
    auto_gen_sampler();
  }
  return old_val;
}

template<class SieveTraits>
void GaussQueue<SieveTraits, false>::auto_gen_sampler()
{
  assert(sampler == nullptr);
  assert(sampler_auto_created == false);
    // Use seed_sampler to init the rng. The other numbers spell "SAMPLER" in ASCII, they are here
    //  because we use the very same seed elsewhere as well (Think of sees_seq as a Hash)

  std::seed_seq seed{83, 65, 77, 80, 76, 69, 82, seed_sampler};
  // For the standalone queue, we do not have a caller sieve and do not know the dimension.
  // Q: Should the initial dimension be set during init from the base rather than at construction?
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing our own sampler. Using GPVSampler")
  sampler = new GPVSampler<SieveTraits, false, std::mt19937_64, std::seed_seq>(seed);
#else
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing our own sampler. Using GPVSamplerExtended")
  sampler = new GPVSamplerExtended<SieveTraits, false, std::mt19937_64, std::seed_seq>(seed, gauss_sieve->get_ambient_dimension() / 3 );
#endif
  sampler_auto_created = true;
//    sampler = new GPVSampler<SieveTraits, false, std::mt19937_64, std::seed_seq>(seed );
//     sampler = new UniformSampler<SieveTraits, false, std::mt19937_64, std::seed_seq>(seed,2);

// Note that the sampler will LATER be initialized from the Sieve itself when the sieve is actually
// run. This late-initialization allows to construct samplers that outlive the main sieve object.
}

/**
  IO:
*/

template<class SieveTraits>
std::ostream &operator<<(std::ostream &os, GaussQueue<SieveTraits,false> const &gauss_queue)
{
  assert(false); // need to take care about internal status of sampler
  return os;
}

template<class SieveTraits>
bool GaussQueue<SieveTraits, false>::serialize_gauss_queue(std::ostream &os) const
{
  assert(false);
  return true;
}

template<class SieveTraits>
bool GaussQueue<SieveTraits, false>::unserialize_gauss_queue(std::istream &is)
{
  assert(false);
  return true;
}


// clang-format on

}  // end namespace GaussSieve

#endif  // include guard
