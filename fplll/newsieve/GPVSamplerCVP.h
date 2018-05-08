// clang-format status: OK

#ifndef GPV_SAMPLER_CVP_H
#define GPV_SAMPLER_CVP_H

#include "DefaultIncludes.h"
#include "LatticeBases.h"
#include "Sampler.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr.h"


namespace GaussSieve
{

template <class SieveTraits, bool MT> class Sieve;

// forward declaration:
template <class SieveTraits, bool MT, class Engine, class Sseq> class GPVSamplerCVP;

template <class SieveTraits, bool MT, class Engine, class Sseq>
class GPVSamplerCVP final : public Sampler<SieveTraits, MT, Engine, Sseq>
{
public:
  using DimensionType = typename SieveTraits::DimensionType;
  using LengthType    = typename SieveTraits::LengthType;
  using RetType       = typename SieveTraits::GaussSampler_ReturnType;

  // clang-format off
  explicit GPVSamplerCVP(Sseq &seq, uint_fast16_t babai, double const _cutoff = 2.0)
      : Sampler<SieveTraits,MT,Engine,Sseq>(seq),
        start_babai(babai),
        cutoff(_cutoff),
        initialized(false),
        static_init_rettype(nullptr),
        static_init_plainpoint(nullptr)
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Constructing GPVSampler.")
  }
  // clang-format on

  virtual SamplerType sampler_type() const override { return SamplerType::GPVCVP_sampler; }
  virtual ~GPVSamplerCVP()
  {
    if (initialized)
    {
      delete static_init_plainpoint;
      delete static_init_rettype;
    }
  }
  virtual inline RetType sample(int const thread = 0) override;

  explicit GPVSamplerCVP(std::istream &is)
    : Sampler<SieveTraits, MT, Engine, Sseq>(is, SamplerType::GPVCVP_sampler),
      //start_babai  // overwritten anyway
      //cutoff(0.0), // overwritten anyway
      initialized(false),
      static_init_rettype(nullptr),
      static_init_plainpoint(nullptr)
  {
      //if (!read_from_stream(is)) throw bad_dumpread("Could not init GPVExtended Sampler from stream");
  }


private:
  inline virtual void custom_init(SieveLatticeBasis<SieveTraits, MT> const &input_basis) override;
 // inline virtual std::ostream &dump_to_stream(std::ostream &os) const override;
 // inline virtual std::istream &read_from_stream(std::istream &is) override;

  // mu_i,j = r_i,j / ||b*_j||^2.. lower triangular matrix
  std::vector<std::vector<double>> mu_matrix;

  std::vector<double> s2pi;           // st. dev. per dimension, squared and divided by pi.
  std::vector<double> maxdeviations;  // [s2pi * cutoff] - for rejection sampling
  DimensionType dim;

  uint_fast16_t lattice_rank;
  uint_fast16_t start_babai;  // start_babai < lattice_rank; use Babai on { b_{start_babai},...b_1}

  double cutoff;
  bool initialized;

  using Sampler<SieveTraits, MT, Engine, Sseq>::sieveptr;
  using Sampler<SieveTraits, MT, Engine, Sseq>::engine;
  std::vector<typename SieveTraits::PlainPoint> basis;

  StaticInitializer<RetType> *static_init_rettype;
  StaticInitializer<typename SieveTraits::PlainPoint> *static_init_plainpoint;
};

}  // namespace GaussSieve

#endif
