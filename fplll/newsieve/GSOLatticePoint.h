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

#define GSOLATTICEPOINT_PARAMS class AbsoluteEntries, class RelativeCoos, class Basis, int rankfixed
#define GSOLATTICEPOINT_FULLNAME GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, rankfixed>

#ifdef DEBUG_SIEVE_CHECKSAMEGSOBASIS
#define CHECKSAMEBASIS(x2) assert(basisptr == x2.basisptr);
#else
#define CHECKSAMEBASIS(x2)
#endif

template<GSOLATTICEPOINT_PARAMS> class GSOLatticePoint;

template <class PresumedBasis> class IsLatticeBasis : public std::false_type {};
template <class Entries, class gEntries, int DimFixed, int RankFixed>
class IsLatticeBasis<LatticeBasis<Entries, gEntries, DimFixed, RankFixed>> : public std::true_type {};

template<GSOLATTICEPOINT_PARAMS>
struct LatticePointTraits<GSOLATTICEPOINT_FULLNAME>
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
#define FOR_FIXED_RANK    template<int rankfixed_copy = rankfixed, TEMPL_RESTRICT_DECL(rankfixed_copy >= 0)>
#define FOR_VARIABLE_RANK template<int rankfixed_copy = rankfixed, TEMPL_RESTRICT_DECL(rankfixed_copy == -1)>
// clang-format on


template<GSOLATTICEPOINT_PARAMS>
class GSOLatticePoint final : public GeneralLatticePoint<GSOLATTICEPOINT_FULLNAME>
{
  static_assert(IsLatticeBasis<Basis>::value == true, "Invalid Basis");
  // constructors are defaulted.
  // Note that copy constructor is automatically deleted, because the parent's is.

  friend StaticInitializer<GSOLATTICEPOINT_FULLNAME>;

public:

  using LatticePointTag = std::true_type;

private:
  using Parent = GeneralLatticePoint<GSOLATTICEPOINT_FULLNAME>;
  // Container type used to store the coefficients of the point wrt an orthonormal basis
  using AbsoluteContainer = ArrayOrVector<AbsoluteEntries, rankfixed>;
  // Coordinates used to store the coeffients of the point wrt the basis
  using CoefficientContainer = ArrayOrVector<RelativeCoos, rankfixed>;
  static MaybeFixed<rankfixed> basis_size;


  // Actual entries:
  AbsoluteContainer onb_coos;  // coordinates wrt to an orthonormal basis
  CoefficientContainer relative_coos;  // coordinates relative to the basis indicated by the static
                                       // pointer
  AbsoluteEntries norm2;  // precomputed norm^2
  std::shared_ptr<Basis> basisptr;

private:

  FOR_FIXED_RANK
  static constexpr MaybeFixed<nfixed> get_rank()
  {
    return MaybeFixed<rankfixed>(rankfixed);
  }

  FOR_VARIABLE_RANK
  static constexpr MaybeFixed<-1> const &get_rank()
  {
    return basis_size;
  }


public:
  // get dimension

  /*
  FOR_FIXED_RANK
  static constexpr MaybeFixed<nfixed> get_dim()
  {
    static_assert(rankfixed_copy == rankfixed, "");  // nfixed_copy defined in FOR_FIXED_DIM
    return MaybeFixed<rankfixed>(rankfixed);
  }

  FOR_VARIABLE_RANK
  static constexpr MaybeFixed<-1> const &get_dim()
  {
    static_assert(rankfixed == -1, "");
    return dim;
  }
  */

  // arithmetic:
  GSOLatticePoint &operator+=(GSOLatticePoint const &x2)
  {
    CHECKSAMEBASIS(x2)
    return Parent::operator+=(*this,x2);
  }
  GSOLatticePoint &operator-=(GSOLatticePoint const &x2)
  {
    CHECKSAMEBASIS(x2)
    return Parent::operator-=(*this, x2);
  }

  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  inline void add_multiply(GSOLatticePoint const &x2, Integer const multiplier)
  {
    CHECKSAMEBASIS(x2)
    static_cast<Parent*>(this)->add_multiply(x2, multiplier);
  }

  template <class Integer, TEMPL_RESTRICT_DECL2(std::is_integral<Integer>)>
  inline void sub_multiply(GSOLatticePoint const &x2, Integer const multiplier)
  {
    CHECKSAMEBASIS(x2)
    static_cast<Parent*>(this)->sub_multiply(x2, multiplier);
  }

  bool operator==(GSOLatticePoint const x2)
  {
    CHECKSAMEBASIS(x2)
    return Parent::operator==(*this, x2);
  }

  FOR_FIXED_RANK
  static constexpr MaybeFixed<nfixed> get_internal_rep_size()
  {
    return MaybeFixed<rankfixed>(rankfixed);
  }

  FOR_VARIABLE_RANK
  static constexpr MaybeFixed<-1> const &get_internal_rep_size()
  {
    return basis_size;
  }

  template <class Arg>
  inline RelativeCoos const &get_internal_rep(Arg &&arg) const { return relative_coos[std::forward<Arg>(arg)]; }
  template <class Arg>
  inline RelativeCoos &get_internal_rep(Arg &&arg) { return relative_coos[std::forward<Arg>(arg)];}
  template <class Arg>
  inline AbsoluteEntries get_absolute_coo(Arg &&arg) const { return onb_coos[std::forward<Arg>(arg)]; }

  auto get_dim() const { return basisptr->get_ambient_dimension(); }

  /*
    Constructor:
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
};

/**
  static member intitalization, compile-time (or startup):
*/
//template<class AbsoluteEntries, class RelativeCoos, class Basis, int nfixed>
//Basis const * GSOLatticePoint<AbsoluteEntries, RelativeCoos, Basis, nfixed>::basisptr = nullptr;
template<GSOLATTICEPOINT_PARAMS>
MaybeFixed<nfixed> GSOLATTICEPOINT_FULLNAME::dim = MaybeFixed<nfixed>(nfixed < 0 ? 0 : nfixed);


/**
  Run-time initalization, via RAII wrapper:
*/

template<GSOLATTICEPOINT_PARAMS>
std::ostream &operator<<(std::ostream &, StaticInitializer<GSOLATTICEPOINT_FULLNAME> const &);
template<GSOLATTICEPOINT_PARAMS>
std::ostream &operator<<(std::istream &, StaticInitializer<GSOLATTICEPOINT_FULLNAME> const &);

template<GSOLATTICEPOINT_PARAMS>
class StaticInitializer<GSOLATTICEPOINT_FULLNAME> final
    : public DefaultStaticInitializer<GSOLATTICEPOINT_FULLNAME>
{
  using Parent = DefaultStaticInitializer<GSOLATTICEPOINT_FULLNAME>;

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
      if (!(new_dim == GSOLATTICEPOINT_FULLNAME::dim))
      {
        throw bad_reinit_static("Trying to reinit static dimension of GSOLatticePoint");
      }
    }
    else
    {
      GSOLATTICEPOINT_FULLNAME::dim = new_dim;
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

  friend std::ostream &operator<<(std::ostream &os, StaticInitializer<GSOLATTICEPOINT_FULLNAME> const &)
  {
    os << "Static Data for GSOLatticePoint ";
    os << ((nfixed < 0) ? "(fixed dim): " : "(variable dim): ");
    os << access_dim();
    return os;
  }
  friend std::istream &operator>>(std::istream &is, StaticInitializer<GSOLATTICEPOINT_FULLNAME> const &init_ob)
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
  static constexpr mystd::add_lvalue_reference_t<decltype(GSOLATTICEPOINT_FULLNAME::dim)> access_dim()
  {
    return GSOLATTICEPOINT_FULLNAME::dim;
  }
};  // end of static initializer



}  // end namespace GaussSieve

#undef FOR_FIXED_DIM
#undef FOR_VARIABLE_DIM
//#undef GSOLATTICEPOINT_PARAMS
//#undef GSOLATTICEPOINT_FULLNAME

#endif  // include guard


