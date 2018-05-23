#warning UNUSED FILE

#ifndef GSO_LATTICE_POINT_H
#define GSO_LATTICE_POINT_H

#include "DefaultIncludes.h"
#include "LatticePointConcept.h"
#include "LatticePointGeneric.h"
#include "LatticeBases.h"

namespace GaussSieve
{

/**
  GSOLatticePoint stores lattice points that are defined wrt a basis.
  Notably, we store a (smart) pointer to a basis of type Basis = GaussSieve::LatticeBasis<...>
  that stores GSO information.
  For each point lp, we store coordinates wrt this basis (i.e. integers relative_coos[i] s.t.
  lp = sum relative_coos[i] * basis[i] ).
  Furthermore, we store (non-integral) coordinates wrt the orthonormal basis defined by the GSO.
  Note: The integral coos are the "primary" defining entries of the point and the non-integral coos
  are considered pre-computed derived data.
  In particular, arithmetic operations are defined with the rel_coo's and we recompute the GSO
  representation (this is done to avoid numeric error growth).

*/

#define GSOLATTICEPOINT_PARAMS \
template <class AbsoluteEntries, class RelativeCoos, class Basis, int nfixed>
#define GSOLATTICEPOINT_WPARAMS \
GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, nfixed>

GSOLATTICEPOINT_PARAMS class GSOLatticePoint;

template <class PresumedBasis> class IsLatticeBasis : public std::false_type {};
template <class Entries, class gEntries, int DimFixed, int RankFixed>
class IsLatticeBasis<LatticeBasis<Entries, gEntries, DimFixed, RankFixed>> : public std::true_type {};

GSOLATTICEPOINT_PARAMS
struct LatticePointTraits<GSOLATTICEPOINT_WPARAMS>
{
  static_assert(IsLatticeBasis<Basis>::value == true, "Invalid Basis");
public:
  using Trait_ScalarProductStorageType      = AbsoluteEntries;
  using Trait_ScalarProductStorageType_Full = AbsoluteEntries;
  using Trait_ExposesCoos                   = std::false_type; // because lp[i] is ambiguous
  //using Trait_CoordinateType                = void;
  //using Trait_Coos_RW                       = std::false_type;
  using Trait_AbsoluteCooType               = AbsoluteEntries;
  using Trait_RepCooType                    = RelativeCoos;
  using Trait_ExposesInternalRep            = std::true_type;
  using Trait_InternalRepLinear             = std::true_type;
  using Trait_InternalRep_RW                = std::true_type;
  using Trait_InternalRepByCoos             = std::false_type;
  using Trait_InternalRepIsAbsolute         = std::false_type;
  //using Trait_AbsoluteCoos                  = std::false_type;
  using Trait_CheapNorm2                    = std::true_type;
  using Trait_CheapNegate                   = std::true_type;
  using Trait_BitApprox                     = std::false_type;
  using Trait_Leveled                       = std::false_type;
};


// local defines, undef at the end of file, used to enable functions only for (non-)fixed dimension.
// clang-format off
#define FOR_FIXED_DIM    template<int nfixed_copy = nfixed, TEMPL_RESTRICT_DECL(nfixed_copy >= 0)>
#define FOR_VARIABLE_DIM template<int nfixed_copy = nfixed, TEMPL_RESTRICT_DECL(nfixed_copy == -1)>
// clang-format on


GSOLATTICEPOINT_PARAMS
class GSOLatticePoint final : public GeneralLatticePoint<GSOLATTICEPOINT_WPARAMS>
{
  static_assert(IsLatticeBasis<Basis>::value == true, "Invalid Basis");
  // constructors are defaulted.
  // Note that copy constructor is automatically deleted, because the parent's is.

  friend StaticInitializer<GSOLATTICEPOINT_WPARAMS>;

public:

  using LatticePointTag = std::true_type;

private:

  // Container type used to store the coefficients of the point wrt an orthonormal basis
  using AbsoluteContainer = ArrayOrVector<AbsoluteEntries, nfixed>;
  // Coordinates used to store the coeffients of the point wrt the basis
  using CoefficientContainer = ArrayOrVector<RelativeCoos, nfixed>;
  static MaybeFixed<nfixed> dim;

public:
  // get dimension
  FOR_FIXED_DIM
  static constexpr MaybeFixed<nfixed> get_dim()
  {
    static_assert(nfixed_copy == nfixed, "");  // nfixed_copy defined in FOR_FIXED_DIM
    return MaybeFixed<nfixed>(nfixed);
  }


  FOR_VARIABLE_DIM
  static constexpr MaybeFixed<-1> const &get_dim()
  {
    static_assert(nfixed == -1, "");
    return dim;
  }

  /*
    Constructors:
  */

  /*

  // clang-format off
  FOR_FIXED_DIM  // This introduces template params, so we do not =default
  constexpr explicit GSOLatticePoint() {}
  // clang-format on

  FOR_FIXED_DIM
  constexpr explicit GSOLatticePoint(MaybeFixed<nfixed> dim) {}

  FOR_VARIABLE_DIM
  explicit GSOLatticePoint(MaybeFixed<nfixed> dim)
      : GSOLatticePoint()
  {
#ifdef DEBUG_SIEVE_LP_MATCHDIM
    assert(dim == get_dim());
#endif
  }

  FOR_VARIABLE_DIM
  explicit GSOLatticePoint()  // make noexcept?
      : onb_coos(get_dim()), relative_coos(get_dim())
  {
    // TODO : assert initialization
  }

  */

  static std::string class_name() { return "Lattice point with coeffs wrt. supplied basis"; }

  AbsoluteEntries get_norm2() const { return norm2; }

  void sanitize() &
  {
    sanitize(compute_sc_product(*this, *this));
  }
  // clang-format on

  // version of the above where we already know norm2. Note that sanitize() just forwards to this.
  void sanitize(AbsoluteEntries const &new_norm2)
  {
    norm2 = new_norm2;
  }

private:
  AbsoluteContainer onb_coos;  // coordinates wrt to an orthonormal basis
  CoefficientContainer relative_coos;  // coordinates relative to the basis indicated by the static
                                       // pointer
  AbsoluteEntries norm2;  // precomputed norm^2
  std::shared_ptr<Basis> basisptr;
};

/**
  static member intitalization, compile-time (or startup):
*/
//template<class AbsoluteEntries, class RelativeCoos, class Basis, int nfixed>
//Basis const * GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, nfixed>::basisptr = nullptr;
GSOLATTICEPOINT_PARAMS
MaybeFixed<nfixed> GSOLATTICEPOINT_WPARAMS::dim = MaybeFixed<nfixed>(nfixed < 0 ? 0 : nfixed);


/**
  Run-time initalization, via RAII wrapper:
*/

GSOLATTICEPOINT_PARAMS
std::ostream &operator<<(std::ostream &, StaticInitializer<GSOLATTICEPOINT_WPARAMS> const &);
GSOLATTICEPOINT_PARAMS
std::ostream &operator<<(std::istream &, StaticInitializer<GSOLATTICEPOINT_WPARAMS> const &);

GSOLATTICEPOINT_PARAMS
class StaticInitializer<GSOLATTICEPOINT_WPARAMS> final
    : public DefaultStaticInitializer<GSOLATTICEPOINT_WPARAMS>
{
  using Parent = DefaultStaticInitializer<GSOLATTICEPOINT_WPARAMS>;

public:
  template <class T, TEMPL_RESTRICT_DECL2(IsArgForStaticInitializer<T>)>
  StaticInitializer(T const &initializer) : StaticInitializer(initializer.dim)
  {
  }

explicit StaticInitializer(MaybeFixed<nfixed> const new_dim)
  {
    assert(Parent::user_count > 0);
    if (Parent::user_count > 1)
    {
      if (!(new_dim == GSOLATTICEPOINT_WPARAMS::dim))
      {
        throw bad_reinit_static("Trying to reinit static dimension of GSOLatticePoint");
      }
    }
    else
    {
      GSOLATTICEPOINT_WPARAMS::dim = new_dim;
    }
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Initializing GSOLatticePoint with nfixed = "
                                      << nfixed << " REALDIM = " << new_dim << " Counter is"
                                      << Parent::user_count)
  }
  ~StaticInitializer()
  {
    DEBUG_SIEVE_TRACEINITIATLIZATIONS("Deinitializing GSOLatticePoint with nfixed = "
                                      << nfixed << " Counter is " << Parent::user_count)
  }

  friend std::ostream &operator<<(std::ostream &os, StaticInitializer<GSOLATTICEPOINT_WPARAMS> const &)
  {
    os << "Static Data for GSOLatticePoint ";
    os << ((nfixed < 0) ? "(fixed dim): " : "(variable dim): ");
    os << access_dim();
    return os;
  }
  friend std::istream &operator>>(std::istream &is, StaticInitializer<GSOLATTICEPOINT_WPARAMS> const &init_ob)
  {
    if (!string_consume(is,"Static Data for GSOLatticePoint")) throw bad_dumpread("Dumpread failure: GSOLatticePoint");
    if (!string_consume(is, (nfixed < 0) ? "(fixed dim):" : "(variable dim):")) throw bad_dumpread("Dumpread failure: GSOLatticePoint");
    mystd::decay_t<decltype(access_dim())> new_dim;
    if (!(is >> new_dim)) throw bad_dumpread("Dumpread failure: GSOLatticePoint");
    if(init_ob.get_user_count() > 1)
    {
      if(new_dim != access_dim()) throw bad_reinit_static("Trying to overwrite static data for GSO lattice point, but more than 1 user");
    }
    access_dim() = new_dim;
    return is;
  }
  explicit StaticInitializer(std::istream &is)
  {
    is >> *this;
  }

  // this allows friends of StaticInitializer to access ExactLatticePoint<ET,nfixed>
  // (since friend is not transitive in C++)
  private:
  static constexpr mystd::add_lvalue_reference_t<decltype(GSOLATTICEPOINT_WPARAMS::dim)> access_dim()
  {
    return GSOLATTICEPOINT_WPARAMS::dim;
  }
};  // end of static initializer



}  // end namespace GaussSieve

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM
#undef GSOLATTICEPOINT_PARAMS
#undef GSOLATTICEPOINT_WPARAMS

#endif  // include guard


