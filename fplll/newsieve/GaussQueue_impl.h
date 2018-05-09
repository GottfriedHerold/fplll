#ifndef GAUSS_QUEUE_CPP
#define GAUSS_QUEUE_CPP

#include "DebugAll.h"
#include "GPVSampler.h"
#include "GPVSamplerExtended.h"
#include "GaussQueue.h"
#include "Sampler.h"
#include "UniformSampler.h"
#include "GPVSamplerCVP.h"

namespace GaussSieve
{

// constructor, single-threaded version
// clang-format off
template <class SieveTraits>
GaussQueue<SieveTraits, false>::GaussQueue(Sieve<SieveTraits, false> *const caller_sieve,
                                          GlobalStaticDataInitializer const &static_data,
                                          int new_seed_sampler,
                                          Sampler<SieveTraits, false> *user_sampler)
    : init_data_type(static_data),
      init_ret_type(static_data),
      main_queue(),
      gauss_sieve(caller_sieve),
      seed_sampler(new_seed_sampler),
      sampler(user_sampler),
      sampler_pointer_owning(false)
// clang-format on
{
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
  assert(caller_sieve == nullptr);
#else
  assert(caller_sieve != nullptr);
#endif
  if (sampler == nullptr)
  {
    auto_gen_sampler(); // create our own sampler and set sampler_pointer_owning
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
  if (sampler_pointer_owning)
  {
    delete sampler;
  }
}

template<class SieveTraits>
auto GaussQueue<SieveTraits, false>::set_sampler(SamplerPtr const &new_sampler_ptr) -> SamplerPtr
{
  if (sampler_pointer_owning)
  {
    delete sampler;
    sampler = nullptr;
    sampler_pointer_owning = false;
  }
  SamplerPtr old_val = sampler;
  sampler = new_sampler_ptr;
  if (sampler == nullptr)
  {
    auto_gen_sampler();
  }
  return old_val;
}

template<class SieveTraits>
void GaussQueue<SieveTraits, false>::auto_gen_sampler()
{
  assert(sampler == nullptr);
  assert(sampler_pointer_owning == false);
    // Use seed_sampler to init the rng. The other numbers spell "SAMPLER" in ASCII, they are here
    //  because we use the very same seed elsewhere as well (Think of sees_seq as a Hash)

  std::seed_seq seed{83, 65, 77, 80, 76, 69, 82, seed_sampler};
  // For the standalone queue, we do not have a caller sieve and do not know the dimension.
  // Q: Should the initial dimension be set during init from the base rather than at construction?
#ifdef DEBUG_SIEVE_STANDALONE_QUEUE
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing our own sampler. Using GPVSampler")
  sampler = new DefaultSamplerTypeStandalone(seed);
#else
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing our own sampler. Using GPVSamplerExtended")
  sampler = new DefaultSamplerType(seed, gauss_sieve->get_ambient_dimension() / 3 );
#endif
  sampler_pointer_owning = true;
//    sampler = new GPVSampler<SieveTraits, false, std::mt19937_64, std::seed_seq>(seed );
//     sampler = new UniformSampler<SieveTraits, false, std::mt19937_64, std::seed_seq>(seed,2);

// Note that the sampler will LATER be initialized from the Sieve itself when the sieve is actually
// run. This late-initialization allows to construct samplers that outlive the main sieve object.
}

/**
  IO: Note that if we use a priority queue, IO just asserts false. The issue is that we cannot
      non-destructively serialize given the current interface.
      Since a regular queue is faster anyway, there is not much point in supporting this.
*/

// NOTE: Partial specialization of non-class function templates is disallowed in C++  :(
// Consequently, this is the code both both MT variants.
// (the friend declaration is fine). If the MT == true / false variants differ, this function
// needs to handoff to a helper class... (or use constexpr if to write two versions inside this
// function)

template<class SieveTraits, bool MT>
std::ostream &operator<< (std::ostream &os, GaussQueue<SieveTraits,MT> const &gauss_queue)
{
  static_assert(MT == false, "This will not work");
#ifndef USE_REGULAR_QUEUE
  assert(false);
#else
  os << "Initialization data: ";
  os << gauss_queue.init_data_type << '\n';
  os << gauss_queue.init_ret_type  << '\n';
  os << "Seed (autogenerated sampler)";
  os << gauss_queue.seed_sampler << '\n';
  os << "Sampler pointer was owning: ";
  os << gauss_queue.sampler_pointer_owning << '\n';
  os << "Internal sampler status: ";
  os << *(gauss_queue.sampler) <<'\n';
  os << "Number of elements in queue: ";
  auto const queue_size = gauss_queue.main_queue.size();
  os << queue_size << '\n';
  os << " [";
  // std::queue does not allow iterator. Consequently, serialization is... difficult.
  // This should work, but is arguably stupid wrt efficiency.
  for (mystd::decay_t<decltype(queue_size)> i = 0; i < queue_size; ++i)
  {
    typename GaussQueue<SieveTraits,false>::DataType buf{std::move(gauss_queue.main_queue.front())};
    gauss_queue.main_queue.pop();
    os << buf << " ";
    gauss_queue.main_queue.push(std::move(buf));
  }
  os << "]";
#endif
  return os;
}


// almost identical to the above, except that we use serialize_lp instead of operator<<
// for the stored lattice point and do slightly different error handling.
template<class SieveTraits>
bool GaussQueue<SieveTraits, false>::serialize_gauss_queue(std::ostream &os) const
{
#ifndef USE_REGULAR_QUEUE
  assert(false);
#else
  os << init_data_type << '\n';
  os << init_ret_type  << '\n';
  os << "Seed";
  os << seed_sampler << '\n';
  os << "Sampler pointer was owning?";
  os << sampler_pointer_owning;
  os << "Sampler type";
  os << static_cast<int>(sampler->sampler_type()) << '\n';
  os << *sampler;
  os << "Number of elements in queue";
  auto const queue_size = main_queue.size();
  os << queue_size;
  os << " [";
  // std::queue does not allow iterator. Consequently, serialization is... difficult.
  // This should work, but is arguably stupid wrt efficiency.
  for (mystd::decay_t<decltype(queue_size)> i = 0; i < queue_size; ++i)
  {
    DataType buf{std::move(main_queue.front())};
    main_queue.pop();
    bool good = buf.serialize_lp(os);
    main_queue.push(std::move(buf));
    if (!good) return false;
  }
  os << "]";
#endif
  return static_cast<bool>(os);
}

template<class SieveTraits>
bool GaussQueue<SieveTraits, false>::unserialize_gauss_queue_main(std::istream &is)
{
#ifndef USE_REGULAR_QUEUE
  assert(false); return false;
#else
  // init_data_type and init_ret_type in unserialize_gauss_queue (without _main)
  if (!string_consume(is,"Seed")) return false;
  if (!(is >> seed_sampler)) return false;
  if (!string_consume(is,"Sampler pointer was owning?")) return false;
  // sampler_pointer_owning == false at this point.
  // Note : We do NOT is >> sampler_pointer_owning, because new ownership of the pointer is not
  // really related to the situation when serializing.
  bool dummy;
  if (!(is >> dummy)) return false;
  if (!string_consume(is,"Sampler type")) return false;
  int tmp;
  if (!(is >> tmp)) return false;
  switch(static_cast<SamplerType>(tmp))
  {
    case SamplerType::GPV_sampler:
    {
      sampler = new GPVSampler<SieveTraits, false, std::mt19937_64, std::seed_seq>(is);
      sampler_pointer_owning = true;
      break;
    }
    case SamplerType::GPVExtended_sampler:
    {
      assert(false); // does not work yet
//      sampler = new GPVSamplerExtended<SieveTraits,false, std::mt19937_64, std::seed_seq>(is);
      sampler_pointer_owning = true;
      break;
    }
    case SamplerType::user_defined:
    {
      if (sampler == nullptr)
      {
        throw std::runtime_error("Cannot initialize unknown Sampler type");
      }
      else
      {
        if(!(is >> *sampler)) return false; // This uses the existing dynamic type of *sampler.
        // This can only happen if we call unserialize when a non-owning pointer that existed
        // beforehand. In this case, we essentially overwrite the object managed by the caller.
        sampler_pointer_owning = false;
      }
      break;
    }
    default:
    {
      throw std::runtime_error("Unknown Sampler type (or not supported yet");
    }
  } // end of case
  if (!string_consume(is,"Number of elements in queue")) return false;
  mystd::decay_t<decltype(main_queue.size())> queue_size;
  if (!(is >> queue_size)) return false;
  if (!string_consume(is,"[")) return false;

  for (decltype(queue_size) i = 0; i < queue_size; ++i)
  {
    DataType buf{};
    if ( !(buf.unserialize_lp(is))) return false;
    main_queue.push(std::move(buf));
  }
  if (!string_consume(is,"]")) return false;
  return true;
#endif
}

template<class SieveTraits>
bool GaussQueue<SieveTraits, false>::unserialize_gauss_queue(std::istream &is)
{
#ifndef USE_REGULAR_QUEUE
  assert(false);
  return false;
#else
  if (sampler_pointer_owning)
  {
    delete sampler;
    sampler = nullptr; // we are going to overwrite it
    sampler_pointer_owning = false;
  }
  if (!(is >> init_data_type)) return false;
  if (!(is >> init_ret_type)) return false;
  return unserialize_gauss_queue_main(is);
#endif
}

template<class SieveTraits>
GaussQueue<SieveTraits, false>::GaussQueue(std::istream &is, Sieve<SieveTraits, false>* caller_sieve)
  : init_data_type(is),
    init_ret_type(is),
    main_queue(),
    gauss_sieve(caller_sieve),
    // seed_sampler(new_seed_sampler),
    sampler(nullptr),
    sampler_pointer_owning(false)
{
#ifndef USE_REGULAR_QUEUE
  assert(false);
#endif
  if (unserialize_gauss_queue_main(is) == false) throw std::runtime_error("Could not create queue from stream");
}


// clang-format on

}  // end namespace GaussSieve

#endif  // include guard
