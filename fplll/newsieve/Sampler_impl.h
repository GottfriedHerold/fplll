// clang-format status: NOT OK (nested ifdefs)

/**
 Implementation for some methods of the sampler. These are in a separate file because they need to
 access the header of the main Sieve class
*/

#ifndef SAMPLER_IMPL_H
#define SAMPLER_IMPL_H

#include "DefaultIncludes.h"
#include "LatticeBases.h"
#include "Sampler.h"
#include "Typedefs.h"

namespace GaussSieve
{

// actually needed, even though destructor is pure virtual as the base class destructor is
// eventually called implicitly.
template <class SieveTraits, bool MT, class Engine, class Sseq>
Sampler<SieveTraits, MT, Engine, Sseq>::~Sampler()
{
}

template <class SieveTraits, bool MT, class Engine, class Sseq>
void Sampler<SieveTraits, MT, Engine, Sseq>::init(
    Sieve<SieveTraits, MT> *const sieve, SieveLatticeBasis<SieveTraits, MT> const &input_basis)
{
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing Sampler:")
  sieveptr = sieve;
#ifndef DEBUG_SIEVE_STANDALONE_SAMPLER
  assert(sieveptr != nullptr);
#else
  assert(sieveptr == nullptr);
#endif

#ifdef DEBUG_SIEVE_STANDALONE_SAMPLER
  engine.init(1);
#else
  engine.init(sieve->get_num_threads());
#endif

#ifdef PROGRESSIVE
#ifndef DEBUG_SIEVE_STANDALONE_SAMPLER
  progressive_rank = sieveptr->get_progressive_rank();
  std::cout << "initialized progressive rank to " << progressive_rank << std::endl;
#else
  // PROGRESSIVE && !STANDALONE :
  progressive_rank = input_basis.lattice_rank;
#endif
#endif

  custom_init(input_basis);
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("Finished Initializing Sampler.")
}

template<class SieveTraits, bool MT, class Engine, class Sseq>
Sampler<SieveTraits, MT, Engine, Sseq>::Sampler(std::istream &is)
  :sieveptr(nullptr) // the other members are overwritten anyway
{
  is >> *this;
}


template <class SieveTraits, bool MT, class Engine, class Sseq>
inline std::ostream &operator<<(std::ostream &os,
                                Sampler<SieveTraits, MT, Engine, Sseq> const &sampler)
{
  os << "Sampler type" << static_cast<int>(sampler.sampler_type()) << '\n';
  os << "Randomness state:" << sampler.engine;
#ifdef PROGRESSIVE
  os << '\n' << "Progressive Rank" << sampler.progressive_rank;
#endif
  return sampler.dump_to_stream(os); // handoff to virtual function.
}

template <class SieveTraits, bool MT, class Engine, class Sseq>
inline std::istream &operator>>(std::istream &is,
                                Sampler<SieveTraits, MT, Engine, Sseq> &sampler)
{

  if (!string_consume(is,"Sampler type")) throw bad_dumpread("Could not read sampler");
  int x;
  if (!(is >> x)) throw bad_dumpread("Could not read sampler");
  if (static_cast<int>(sampler.sampler_type()) != x) throw bad_dumpread("Wrong type of sampler");
  if (!string_consume(is,"Randomness state:")) throw bad_dumpread("Could not read sampler");
  if (!(is >> sampler.engine)) throw bad_dumpread("Could not read sampler");
#ifdef PROGRESSIVE
  if (!string_consume(is,"Progressive Rank")) throw bad_dumpread("Could not read sampler");
  if (!(is >> sampler.progressive_rank)) throw bad_dumpread("Could not read sampler");
#endif
  return sampler.read_from_stream(is);
}

// explicit instantiation of template
// template class MTPRNG<std::mt19937_64, false, std::seed_seq>;

}  // end namespace GaussSieve

#endif  // include guards
