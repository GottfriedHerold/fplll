/**
This file provides an interface for a PRNG in a possibly multi-threaded setting.
Notably, we use a master seed to initialize #threads many child RNGs that operate independently.
This is done to make sure that (at least if we don't change basis), the set of vectors that are
sampled is at least somewhat consistent among different runs.

This file also contains some basic rejection sampling routines for sampling Gaussians on Z.
*/

#ifndef MTPRNG_H
#define MTPRNG_H

#include "DefaultIncludes.h"
#include "SieveUtility.h"
#include "Typedefs.h"

namespace GaussSieve
{
// forward declaration:
// wrapper around (a vector of) random number engines of type Engine
// This is used to unify the single and multi-threaded case
template <class Engine, bool MT, class Sseq> class MTPRNG;
template <class Engine, bool MT, class Sseq>
std::ostream &operator<<(std::ostream &, MTPRNG<Engine, MT, Sseq> const &);
template <class Engine, bool MT, class Sseq>
std::istream &operator>>(std::istream &, MTPRNG<Engine, MT, Sseq> &);

/**

 These functions sample from a discrete Gaussian distribution with parameter s and center c
 on (the 1-dimensional lattice given by) the integers Z.

 We cutoff the Gaussian at s*cutoff. This means that the distribution is discrete on a subset of
 Z with output probability for x being proportional to exp(-pi(x-c)^2/s^2). Note the scaling by pi
 in the exponent.

 For reasons of numerical stability, center should not be very large in absolute value (it is
 possible to reduce to |center|<1 anyway), such that  center +/- cutoff * s does not overflow.
 Z must be an integral POD type (e.g. short, int, long, long long).
 center needs to be representable as an (exact) double.
 We do NOT support mpz_t here! The output will take the role of coefficients wrt a given basis.

 We only support double for the floating point numbers. For sieving algorithms, there is no really
 good reason for now to support different precisions, as sampling does not dominate anyway.
 Furthermore, the algorithm is not very sensitive to the quality of the samples.

 Note: if one ever wants to have higher precision, one also needs to adjust the PRNGs to actually
 output high precision. (std::... only supports double)

 engine is supposed to be a random number engine (as defined by the STL).

 The variant sample_z_gaussian_VMD takes
 s2pi = s^2 / pi and maxdeviation = cutoff * s as parameters.

 The implementation is just a simple rejection sampling.
 Note that for the sieve, we do NOT need high quality Gaussians.
*/

// the template argument Z (e.g. int) has to be explicitly provided.
template <class Z, class Engine>
Z sample_z_gaussian(double s, double const center, Engine &engine, double const cutoff);

template <class Z, class Engine>
Z sample_z_gaussian_VMD(double const s2pi, double const center, Engine &engine,
                        double const maxdeviation);

// Samples uniformly between 0 and max_val (inclusive)
template <class Engine> unsigned int sample_uniform(unsigned int max_val, Engine &engine);


/**
The class MTPRNG is just a wrapper around a PRNG engine to facilitate switching to multi-threaded.
Due to the fact that we support multi-threading, MTPRNG<Engine,true,.> is a wrapper around
a vector of Engines, whereas, MTPRNG<Engine,false,.> is a wrapper around Engine.
reseed seeds *all* engines. Use rnd(thread-id) to obtain the underlying Engine.
Thread-safety: Init and reseed are not thread-safe. Concurrent calls to rnd are fine.
You may concurrently call rnd and use init to increase (but not decrease) the number of threads.

Randomness path:
The global (master) seed that is input to the constructor resp. to reseed is used to create
20x32 = 640 bit per-thread-seeds for each actual Engine.

The output from theses engine(s) is then accessed using rnd(thread-number) with
0<=thread-number < num_threads.
Note that rnd(thread-number) returns a reference to the engine. To get actual random data,
use rnd(thread-number)() or feed the engine to a distribution, e.g.
std::uniform_int_distribution<int> six_sided_die(1,6)
int result = six_sided_die(rnd(thread-id));

The usage syntax is the same for the single-threaded case to ensure consistency.
(in particular, rnd takes a thread-id, which should be zero). The randomness path is the same,
so for a given master seed, rnd(0) does not depend on whether the single- or multi-threaded variant
is used.
Note that for obtaining the per-thread seeds from the master seeds, we use a fixed Mersenne twister
engine and not the engine given as template parameter.
*/



/************************
 * Multi-threaded case. *
 ************************/

template <class Engine, class Sseq> class MTPRNG<Engine, true, Sseq>
{
public:
  // clang-format off
  // constructs an uninitialized MTPRNG
  // note that _seq is actually changed and unusable after the call!
  explicit MTPRNG(Sseq &_seq);

  // clang-format on
  inline void reseed(Sseq &_seq);

  /**
  will make sure at least _num_threads engines are actually running, starting new ones as desired.
  Will never reseed/restart already running engines. Using init to reduct the number of threads
  "suspends" the corrensponding rngs (calling them is disallowed).
  Reducing the number of threads and increasing it back saves the random state (unless we reseed).
  */
  inline void init(unsigned int const _num_threads);

  inline Engine &rnd(unsigned int const thread)
  {
#ifdef DEBUG_SIEVE_MTPRNG_THREAD_RANGE
    assert(thread < num_threads);  // might actually work if num_threads < engines.size()
#endif
    return engines[thread];
  }

// operators<< and >> just relay to those:

  bool serialize_mtprng(std::ostream &os) const;
  bool unserialize_mtprng(std::istream &is);
  explicit MTPRNG(std::istream &is);
  bool operator==(MTPRNG const &other) const;
  bool operator!=(MTPRNG const &other) const;

private:
  // seeded with initial seq and consecutively used to seed the children PRNGs.
  std::mt19937_64 seeder; // we use this exactly engines.size() many times after (re-)seeding,
                          // irrespective of changes to num_threads
  std::vector<Engine> engines;
  // number of "active" engines. May differ from size of the vector.
  // In particular, num_threads = 0 means uninitialized.
  unsigned int num_threads;
  /**
     number of 32bit values to use as seeds for the underlying engine(s). Technically, we could use
     state_size if the engine provides it, but not even all default engines do.
  */
  static unsigned int constexpr seed_length = 20;
};

/************************
 * Single-threaded case *
 ************************/

// just wrapper around Engine
template <class Engine, class Sseq> class MTPRNG<Engine, false, Sseq>
{
public:
  explicit MTPRNG(Sseq &_seq) : engine()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing ST RNG Engine.")
    reseed(_seq);
  }
  inline void reseed(Sseq &_seq);

  // does nothing. The (ignored) argument corresponds is to ensure a consistent interface to
  // to the multithreaded case, where it is the number of threads.
  FORCE_INLINE inline void init(unsigned int const = 1)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing Single-Threaded RNG Engines.")
  }

  bool serialize_mtprng(std::ostream &os) const
  {
    os << engine;
    return static_cast<bool>(os);
  }

  bool unserialize_mtprng(std::istream &is)
  {
    is >> engine;
    return static_cast<bool>(is);
  }

  explicit MTPRNG(std::istream &is)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing ST RNG from stream")
    if (!(unserialize_mtprng(is))) throw bad_dumpread("Could not deserialize RNG");
  }

  // Argument is number of thread. It is ignored.
  // Note: non-const constexpr (in C++14) is exactly what we need.
  CPP14CONSTEXPR inline Engine &rnd(unsigned int const = 0) { return engine; }
  bool operator==(MTPRNG const &other) const { return engine == other.engine; }
  bool operator!=(MTPRNG const &other) const { return engine != other.engine; }

private:
  Engine engine;
  /**
  number of 32bit values to use as seeds for the underlying engine(s). Technically, we could
  use state_size if the engine provides it, but not even all default engines do.
    */
  static unsigned int constexpr seed_length = 20;
};  // End of MTPRNG


}  // end of namespace GaussSieve

#include "MTPRNG_impl.h"

#endif
