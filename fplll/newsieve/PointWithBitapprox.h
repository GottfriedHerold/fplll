#ifndef POINT_WITH_BITAPPROX_H
#define POINT_WITH_BITAPPROX_H

#include "DefaultIncludes.h"
#include "GlobalBitApproxData.h"


#define T_IS_ELP static_assert(std::is_same<T,ELP>::value,"Wrong template argument")

// Forward traits:
template<class ELP, class CooSelection>
class LatticePointTraits< AddBitApproximationToLP<ELP,CooSelection> >
{
static_assert(IsALatticePoint<ELP>::value,"ELP is no lattice point");
public:
// forwarding traits from ELP, except for BitApprox
  using Trait_ScalarProductStorageType = Get_ScalarProductStorageType<ELP>;
  using Trait_ScalarProductStorageType_Full  = MakeLeveledScalar<Get_ScalarProductStorageType<ELP>>;
  using Trait_CoordinateType          = Get_CoordinateType<ELP>;
  using Trait_AbsoluteCoos            = Get_AbsoluteCooType<ELP>;
  using Trait_RepCooType              = Get_RepCooType<ELP>;
  using Trait_ExposesCoos             = NormalizeTrait<Has_ExposesCoos<ELP>>;
  using Trait_Coos_RW                 = NormalizeTrait<Has_Coos_RW<ELP>>;
  using Trait_ExposesInternalRep      = NormalizeTrait<Has_ExposesInternalRep<ELP>>;
  using Trait_InternalRepLinear       = NormalizeTrait<Has_InternalRepLinear<ELP>>;
  using Trait_InternalRep_RW          = NormalizeTrait<Has_InternalRep_RW<ELP>>;
  using Trait_InternalRepByCoos       = NormalizeTrait<Has_InternalRepByCoos<ELP>>;
  using Trait_InternalRepIsAbsolute   = NormalizeTrait<Has_InternalRepIsAbsolute<ELP>>;
  using Trait_CheapNorm2              = NormalizeTrait<Has_CheapNorm2<ELP>>;
  using Trait_CheapNegate             = NormalizeTrait<Has_CheapNegate<ELP>>;
  using Trait_BitApprox               = std::true_type;
  using Trait_Leveled                 = NormalizeTrait<Has_Leveled<ELP>>;
  using Trait_ApproxLevel             = Get_ApproxLevel<ELP>;
  // This would not work as expected:
  static_assert(Has_BitApprox<ELP>::value == false,"Trying to add 2 bitapproximations");
};

template<class ELP, class CooSelection>
class AddBitApproximationToLP
: public GeneralLatticePoint<AddBitApproximationToLP<ELP,CooSelection>>
{
  public:
  using LatticePointTag         = std::true_type;
  using Myself = AddBitApproximationToLP<ELP,CooSelection>;
  using SimHashBlock = typename CooSelection::SimHashBlock;
  using SimHashes    = typename CooSelection::SimHashes;

  using PlainCooType  = Get_CoordinateType<ELP>;  // may be void
  using RepCooType    = Get_RepCooType<ELP>;
  using ScalarProductStorageType = Get_ScalarProductStorageType<ELP>;
  using ScalarProductStorageType_Full = Get_ScalarProductStorageType_Full<Myself>;

  private:
  ELP elp;
  SimHashes sim_hashes;
  public:
  constexpr       AddBitApproximationToLP(ELP const &v) = delete;
  CPP14CONSTEXPR  AddBitApproximationToLP(ELP      &&v)
      : elp(std::move(v)),
        sim_hashes(GlobalBitApproxData<CooSelection>::coo_selection.compute_all_bitapproximations(*this)) {}

  constexpr       operator ELP() const & {return elp;}
  CPP14CONSTEXPR  operator ELP()      && {return std::move(elp);}

  void update_bit_approx()
  {
    sim_hashes = GlobalBitApproxData<CooSelection>::coo_selection.compute_all_bitapproximations(*this);
  }

  SimHashBlock const & access_bitapproximation(unsigned int level) const
  {
    return sim_hashes[level];
  }


// forward functionality of ELP. Note that the this is *not* captured by the conversion operators,
// because MakeLeveledLatticePoint<ELP> is derived from GeneralLatticePoint, and those defaults will have
// precendence over any potential conversions.

  static std::string class_name() { return ELP::class_name() + " with Bitapproximations"; }

  template<class T=ELP, class Arg>
  inline PlainCooType &      operator[](Arg &&arg)
  {
    T_IS_ELP; static_assert(Has_ExposesCoos<T>::value && Has_Coos_RW<T>::value,"");
    return elp[std::forward<Arg>(arg)];
  }
  template<class T=ELP, class Arg>
  inline PlainCooType const& operator[](Arg &&arg) const
  {
    T_IS_ELP; static_assert(Has_ExposesCoos<T>::value,"");
    return elp[std::forward<Arg>(arg)];
  }

// operators<,>,<=, >= : No overloads. Defaults is correct.
// forward +=,*=,-=,unary-
  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  inline Myself& operator+=(LatP2 &&x2) { elp+=std::forward<LatP2>(x2); return *this; }
  template<class LatP2, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP2>)>
  inline Myself& operator-=(LatP2 &&x2) { elp-=std::forward<LatP2>(x2); return *this; }
  template<class Multiplier>
  inline Myself& operator*=(Multiplier &&x2) { elp*=std::forward<Multiplier>(x2); return *this; }

  inline Myself operator-() &&   { return static_cast<Myself>(-std::move(elp)); }
  inline bool operator==(Myself const &x2) const { return elp == x2.elp; }
  inline bool operator==(ELP const &x2) const { return elp=x2; }

  // forward get_internal_rep_size, get_internal_rep
  template<class T=ELP>
  CPP14CONSTEXPR inline auto get_internal_rep_size() const -> decltype( std::declval<T>().get_internal_rep_size() )
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value,"");
    return elp.get_internal_rep_size();
  }
  template<class T=ELP, class Arg>
  inline RepCooType const & get_internal_rep(Arg &&arg) const
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value,"");
    return elp.get_internal_rep(std::forward<Arg>(arg));
  }
  template<class T=ELP, class Arg>
  inline RepCooType & get_internal_rep(Arg &&arg)
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value && Has_InternalRep_RW<T>::value,"");
    return elp.get_internal_rep(std::forward<Arg>(arg));
  }

  // forward get_absolute_coo
  template<class Arg> auto inline get_absolute_coo(Arg &&arg) const
    -> decltype( std::declval<ELP>().get_absolute_coo( std::declval<Arg &&>() ))
  { return elp.get_absolute_coo(std::forward<Arg>(arg));  }

  // forward get_dim
  CPP14CONSTEXPR auto inline get_dim() const -> decltype( std::declval<ELP>().get_dim() ) { return elp.get_dim(); }

  // forward write_lp_to_stream
  inline std::ostream& write_lp_to_stream(std::ostream &os, bool const include_norm2=true, bool const include_approx =true) const
  {
    elp.write_lp_to_stream(os, include_norm2,include_approx);
    if (include_approx)
    {
      os << sim_hashes;
    }
    return os;
  }

  // forward write_lp_rep_to_stream
  template<class T=ELP>
  inline std::ostream& write_lp_rep_to_stream(std::ostream &os) const
  {
    T_IS_ELP; static_assert(Has_ExposesInternalRep<T>::value,"");
    return elp.write_lp_rep_to_stream(os);
  }

//
//  //TODO: read_from_stream
//

  void fill_with_zero() { elp.fill_with_zero(); }
  void make_negative()  { elp.make_negative(); }
  bool is_zero() const { return elp.is_zero(); }

  AddBitApproximationToLP make_copy() const & { return static_cast<AddBitApproximationToLP>(elp.make_copy()); }

// TODO: More overloads?
  void sanitize() { elp.sanitize(); }
  template<class Arg>
  void sanitize( Arg &&arg) { elp.sanitize(std::forward<Arg>(arg)); }

// forward get_norm2()
  inline auto get_norm2() const -> decltype (std::declval<ELP>().get_norm2())
  { return elp.get_norm2(); }

  template<unsigned int level>
  inline auto get_norm2_at_level() const -> decltype (std::declval<ELP>().get_norm2_at_level<level>() )
  { return elp.get_norm2_at_level<level>(); }

  inline auto get_norm2_full() const -> decltype (std::declval<ELP>().get_norm2_full())
  { return elp.get_norm2_full(); }

  // forward scalar products:
  template<class Arg>
  inline auto do_compute_sc_product(Arg &&arg) const
    -> decltype( std::declval<ELP>().do_compute_sc_product(std::declval<Arg>() ) )
  {
    return elp.do_compute_sc_product(std::forward<Arg>(arg));
  }

  template<unsigned int level, class Arg>
  inline auto do_compute_sc_product_at_level(Arg &&arg) const
    -> decltype( std::declval<ELP>().do_compute_sc_product_at_level<level>(std::declval<Arg>() ) )
  {
    return elp.do_compute_sc_product_at_level<level>(std::forward<Arg>(arg));
  }

  template<class Arg>
  inline auto do_compute_sc_product_full(Arg &&arg) const
    -> decltype( std::declval<ELP>().do_compute_sc_product_full(std::declval<Arg>() ) )
  {
    return elp.do_compute_sc_product_full(std::forward<Arg>(arg));
  }
};

// static initializer
template<class ELP, class CooSelection>
class StaticInitializer<AddBitApproximationToLP<ELP,CooSelection>>
final : public DefaultStaticInitializer<AddBitApproximationToLP<ELP,CooSelection>>
{
  StaticInitializer<ELP>  const init_elp;
  public:
  template<class X>
  explicit StaticInitializer(X &&init_arg) : init_elp(std::forward<X>(init_arg))
  {
    static_assert(IsArgForStaticInitializer<mystd::decay_t<X>>::value,"");
  }
};

#undef T_IS_ELP

#endif // POINT_WITH_BITAPPROX_H