#ifndef MTPRNG_IMPL_H
#define MTPRNG_IMPL_H
#ifndef MTPRNG_H
#error this file should be included from MTPRNG.h
#endif

namespace GaussSieve
{

template<class Engine, bool MT, class Sseq>
std::ostream &operator<<(std::ostream &os, MTPRNG<Engine, MT, Sseq> const &mtprng)
{
  if (mtprng.serialize_mtprng(os) == false) throw std::runtime_error("Could write MTPRNG");
  return os;
}

template<class Engine, bool MT, class Sseq>
std::istream &operator>>(std::istream &is, MTPRNG<Engine, MT, Sseq> &mtprng)
{
  if (mtprng.unserialize_mtprng(is) == false) throw bad_dumpread("Could not read in MTPRNG");
  return is;
}

template<class Engine, class Sseq>
MTPRNG<Engine, true, Sseq>::MTPRNG(Sseq &_seq) : engines(0), num_threads(0)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing (yet-uninitialized) MT RNG engines.")
    reseed(_seq);
  }

template<class Engine, class Sseq>
bool MTPRNG<Engine, true, Sseq>::serialize_mtprng(std::ostream &os) const
{
  os << "MTPRNG ";
  os << seeder << " ";
  os << "Number of Engines" << engines.size() << '\n' << "[";
  for(typename std::vector<Engine>::size_type i = 0; i < engines.size(); ++i)
  {
    os << engines[i] << " , ";
  }
  os << "]" << '\n';
  os << "threads active " << num_threads;
  return static_cast<bool>(os);
}

template<class Engine, class Sseq>
bool MTPRNG<Engine, true, Sseq>::unserialize_mtprng(std::istream &is)
{
  if (!string_consume(is,"MTPRNG")) return false;
  if (!(is >> seeder)) return false;
  if (!string_consume(is,"Number of Engines")) return false;
  engines.clear();
  assert(engines.size()==0);
  typename decltype(engines)::size_type num;
  if (!(is >> num)) return false;
  if (!string_consume(is,"[")) return false;
  for(decltype(num) i = 0; i < num; ++i)
  {
    Engine eng;
    if (!(is >> eng)) return false;
    if (!string_consume(is, ",")) return false;
    engines.push_back(std::move(eng));
  }
  if (!string_consume(is,"]")) return false;
  if (!string_consume(is,"threads active")) return false;
  if (!(is >> num_threads)) return false;
  return true;
}

template<class Engine, class Sseq>
MTPRNG<Engine, true, Sseq>::MTPRNG(std::istream &is)
{
  if (unserialize_mtprng(is) == false) throw bad_dumpread("Failed to construct MTPRNG");
}

template<class Engine, class Sseq>
bool MTPRNG<Engine, true, Sseq>::operator==(MTPRNG const &other) const
{
  return (num_threads == other.num_threads) && (seeder == other.seeder)
      && (engines == other.engines); // std::vector compares size && element-wise
}

template<class Engine, class Sseq>
bool MTPRNG<Engine, true, Sseq>::operator!=(MTPRNG const &other) const
{
  return !((*this) == other);
}

// clang-format off
template <class Engine, class Sseq>
inline void MTPRNG<Engine, true, Sseq>::reseed(Sseq &_seq)
// clang-format on
{
  seeder.seed(_seq); // seed the pre-RNG seeder whose output is used to seed the engines.
  unsigned int old_threads = num_threads;
  num_threads              = 0;
  engines.clear();
  assert(engines.size()==0);
  init(old_threads);  // will restart all engines.
}

template <class Engine, class Sseq>
inline void MTPRNG<Engine, true, Sseq>::init(unsigned int const new_num_threads)
{
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("(Re-)initializing Multithreaded RNG Engines")
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("setting number of threads as" << new_num_threads)

  num_threads = new_num_threads;
  if (new_num_threads <= engines.size())
  {
    return;
  }
  unsigned int const old_size = engines.size();
  engines.resize(num_threads);
  engines.shrink_to_fit();
  uint32_t per_engine_seed[seed_length];
  // else initialize remaining threads
  for (unsigned int i = old_size; i < num_threads; ++i)
  {

    for (unsigned int j = 0; j < seed_length; ++j)
    {
      per_engine_seed[j] = seeder();
    }
    std::seed_seq per_engine_see_seq(per_engine_seed, per_engine_seed + seed_length);
    engines[i].seed(per_engine_see_seq);
  }
}

/**
  Single-theaded implementation details
*/

// clang-format off
template <class Engine, class Sseq>
inline void MTPRNG<Engine, false, Sseq>::reseed(Sseq &_seq)
// clang-format on
{
  // preprocess _seq in a way that matches the multithreaded case
  std::mt19937_64 seeder(_seq); // not std::move here (
  uint32_t per_engine_seed[seed_length];
  for (unsigned int j = 0; j < seed_length; ++j)
  {
    per_engine_seed[j] = seeder();
  }
  std::seed_seq derived_seed_seq(per_engine_seed, per_engine_seed + seed_length);
  engine.seed(derived_seed_seq);
}


/**
  implementation of the auxiliary samplers, not part of MTPRNG
*/

template <class Z, class Engine>
inline Z sample_z_gaussian(double s, double const center, Engine &engine, double const cutoff)
{
  static_assert(std::is_integral<Z>::value,
                "Return type for sample_z_gaussian must be POD integral type.");

  // maximum deviation of the Gaussian from the center. Note that maxdev may be 1.
  Z const maxdev = static_cast<Z>(std::ceil(s * cutoff));

  // uniform distribution on the set of possible outputs.
  std::uniform_int_distribution<Z> uniform_in_range(std::floor(center - maxdev),
                                                    std::ceil(center + maxdev));

  // defaults to value from [0,1), used in rejection sampling to make the decision.
  std::uniform_real_distribution<double> rejection_test(0.0, 1.0);

  // closest int to center, i.e. most likely value.
  Z const closest_int = std::round(center);

  // negative squared distance to most likely value.
  // This is used to scale up the Gaussian weight function s.t. it is 1 at the most likely value.
  double const adj = -(center - closest_int) * (center - closest_int);

  s = s * s / pi;  // overwriting (local copy of) s. We only care about s^2/pi anyway.


  // use rejection sampling
  while (true)
  {
    Z result    = uniform_in_range(engine);  // sample uniform result.
    double dist = result - center;
    // compute Gaussian weight. std::fma(dist,dist,adj) computes dist^2 + adj = (result-center)^2  -
    // MIN{(result-center)^2 | result integral}.
    //(fma = fused-multiply-add)
    // Recall that s was overwritten to be s^2/pi.

    // Note that the argument of the exp-function might be a tiny positive value due to numeric
    // error
    //(Even if result==closest_int, adj = ROUND((closest_int-center)^2), the computation of
    // std::fma(dist,dist,adj) does not round the intermediate dist^2, leading to a non-zero
    // argument)
    // In particular, it is conceivable that floating point underruns occur in the std::fma - call.
    // Furthermore, if cutoff is large or if s<<1 (in this case, the issue is the rounding when we
    // determined the range), the argument to exp can be extremely small, leading to further
    // potential underruns.
    // We do not care about this for now...

    if (rejection_test(engine) < std::exp(-std::fma(dist, dist, adj) / s))
    {
      return result;
    }
  }
}

// Version taking in s^2/pi and (absolute) maximum deviation. Works just as above.
template <class Z, class Engine>
inline Z sample_z_gaussian_VMD(double const s2pi, double const center, Engine &engine,
                               double const maxdeviation)
{
  static_assert(std::is_integral<Z>::value,
                "Return type for sample_z_gaussian must be POD integral type.");
  std::uniform_int_distribution<Z> uniform_in_range(std::floor(center - maxdeviation),
                                                    std::ceil(center + maxdeviation));
  std::uniform_real_distribution<double> rejection_test(0.0, 1.0);
  Z const closest_int = std::round(center);
  double const adj    = -(center - closest_int) * (center - closest_int);


  while (true)
  {
    Z result    = uniform_in_range(engine);  // sample uniform result.
    double dist = result - center;

    if (rejection_test(engine) < std::exp(-std::fma(dist, dist, adj) / s2pi))
    {
      return result;
    }
  }
}

// Samples uniformly at random from the interval [0, max_val]
// clang-format off
template <class Engine>
inline unsigned int sample_uniform(unsigned int max_val, Engine &engine)
// clang-format on
{
  std::uniform_int_distribution<int> uniform_in_range(0, max_val);
  return uniform_in_range(engine);
}


// Samples non-integer Gaussians, no rejection sampling

template<class FP, class Engine>
inline FP sample_fp_gaussian(double const s2pi, double const mean, Engine &engine)
{
  // std::normal_distribution is restricted to these types.
  static_assert(std::is_floating_point<FP>::value, "Only use for floats");
  // FP must not be a reference or const.
  static_assert(std::is_same<FP,mystd::decay_t<FP>>::value, "");

  std::normal_distribution<FP> gaussian_distribution(mean, s2pi);
  FP result = gaussian_distribution(engine);
  return result;

}



} // end namespace GaussSieve

#endif
