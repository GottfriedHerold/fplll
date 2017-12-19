#ifndef GAUSS_VECTOR_BITAPPROX_H
#define GAUSS_VECTOR_BITAPPROX_H

#include "DefaultIncludes.h"
#include "LatticePointConcept.h"
#include "SieveUtility.h"
#include "SimHash.h"
#include "Typedefs.h"

namespace GaussSieve
{

template <class SieveTraits, bool MT> class GaussVectorWithBitApprox;
template <class SieveTraits, bool MT> class GaussVectorIterator;



// we directly store an approximation (which is actually exact in most cases) to norm2 in the list
// nodes. The class storing this should have a fixed size in memory (no dynamic allocation!)
// and NOT depend on the input data types.

using SimHashApproxNorm2 = double;

namespace GaussListDetails
{
// clang-format off
template <class SieveTraits>
class STVecNode
{
  friend GaussVectorWithBitApprox<SieveTraits, false>;
  friend GaussVectorIterator<SieveTraits, false>;

  // retrieve typedefs to avoid having to write long names.
private:  // shorthands
  using CoordinateSelectionUsed = typename SieveTraits::CoordinateSelectionUsed;
  using SimHashes               = typename CoordinateSelectionUsed::SimHashes;
  using GlobalSimHashClass      = GlobalBitApproxData<CoordinateSelectionUsed>;
  using MainStoragePointer      = typename SieveTraits::MainStoragePointer;
  // This would store 2 independently maintained sim_hashes.
  static_assert(Has_BitApprox<typename SieveTraits::MainStorageType>::value == false, "");

public:
  // we should never need to copy.
  STVecNode(STVecNode const  &)          = delete;
  STVecNode(STVecNode       &&) noexcept = default;
  // clang-format on

  // We construct a STNode from a lattice point.
  // We differentiate whether the input lattice point already stores a sim_hash
  // (if no, we compute a sim_hash, if yes, we take it from the point).
  // Note that the initializer list also (tries to) convert the point, so if the argument
  // type differs from GaussList_StoredPoint, it needs to be convertible
  // (AddBitApproximationToPoint (cf. PointWithApproximation.h) ensures the latter)
  // The !is_reference<Arg> in the restrictions makes this template only valid for actual rvalues.

  // Variant for Arg's without SimHashes
  template <class Arg, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<Arg>>,
                                            mystd::negation<Has_BitApprox<mystd::decay_t<Arg>>>)>
  explicit STVecNode(Arg &&arg) noexcept  // Note: new might actually throw.
      : bit_approximations(GlobalSimHashClass::coo_selection.compute_all_bitapproximations(arg)),
        approx_norm2(convert_to_double(arg.get_norm2())), //ptr_to_exact(nullptr)
        ptr_to_exact(new typename SieveTraits::MainStorageType(std::move(arg)))
  {
    static_assert(std::is_reference<Arg>::value == false, "Use move semantics");
  }

  // Variant for Arg's with SimHashes
  template <class Arg, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<Arg>>,
                                            Has_BitApprox<mystd::decay_t<Arg>>)>
  explicit STVecNode(Arg &&arg) noexcept  // Note: new might actually throw...
      : bit_approximations(std::move(arg).take_bitapproximations()),
        approx_norm2(convert_to_double(arg.get_norm2())),
        ptr_to_exact(new typename SieveTraits::MainStorageType(std::move(arg)))
  {
    static_assert(std::is_reference<Arg>::value == false, "Use move semantics");
  }

  // currently unused: Takes a pointer to a Lattice Point (of the correct type)
  //                   This version has the advantage of not reallocating memory.
  //                   NOTE: THIS TAKES OWNERSHIP OF THE POINTER
  explicit constexpr STVecNode(typename SieveTraits::MainStoragePointer &&point_ptr) noexcept
      : bit_approximations(
            GlobalSimHashClass::coo_selection.compute_all_bitapproximations(*point_ptr)),
        approx_norm2(convert_to_double(point_ptr->get_norm2())),
        ptr_to_exact(point_ptr)
  {
  }

  ~STVecNode() noexcept { delete ptr_to_exact; }

  void update_bitapproximations()
  {
    bit_approximations = GlobalSimHashClass::coo_selection.compute_all_bitapproximations(*ptr_to_exact);
  }

private:
  SimHashes bit_approximations;
  SimHashApproxNorm2 approx_norm2;  // currently a double, consider making it a float
  MainStoragePointer ptr_to_exact;  // owning pointer???

public:
  // we may keep the points sorted according to length. We define operator< to be able to use use
  // sort() from std::
  // (Note that comparions for the lattice points *ptr_to_exact are by norm2)

  // TODO: Consier approx_norm2 instead???
  bool operator<(STVecNode const &other) const { return *ptr_to_exact < *(other.ptr_to_exact); }
};


}  // end namespace Details






template<class SieveTraits>
class GaussVectorWithBitApprox<SieveTraits, false>
{
private:
//  using StoredPoint             = typename SieveTraits::GaussList_StoredPoint;
//  using ReturnType              = typename SieveTraits::GaussList_ReturnType;
  using CoordinateSelectionUsed = typename SieveTraits::CoordinateSelectionUsed;
  using MainStorageType         = typename SieveTraits::MainStorageType;

  using SimHashBlock            = typename CoordinateSelectionUsed::SimHashBlock;
  using SimHashes               = typename CoordinateSelectionUsed::SimHashes;

  using Iterator                = GaussVectorIterator<SieveTraits, false>;  // custom iterator class below
  // clang-format on

  friend Iterator;
  // The class is essentially just a wrapper around UnderlyingContainer
  using UnderlyingContainer         = std::vector<GaussListDetails::STVecNode<SieveTraits>>;
  using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;

public:
  // clang-format off
  GaussVectorWithBitApprox()                                            = delete;
  GaussVectorWithBitApprox(GaussVectorWithBitApprox const &)            = delete;
  GaussVectorWithBitApprox(GaussVectorWithBitApprox &&)                 = delete;
  GaussVectorWithBitApprox &operator=(GaussVectorWithBitApprox const &) = delete;
  GaussVectorWithBitApprox &operator=(GaussVectorWithBitApprox &&)      = delete;
  // clang-format on

  // Sole constructor. The argument is used to initialize the static data of the used lattice point
  // classes.
  // NOTE / TODO: Making this noexcept means we have to catch reinitializations in the caller.
  explicit GaussVectorWithBitApprox(GlobalStaticDataInitializer const &static_data) noexcept
      : init_stored_point_type(static_data),
        actual_vector()
  {
  }

  // behaves like cbegin, cend from STL containers, i.e. gives const-iterator to begin/end.
  // clang-format off
  CPP14CONSTEXPR Iterator begin() noexcept { return actual_vector.begin(); }
  CPP14CONSTEXPR Iterator end()   noexcept { return Iterator{actual_vector.end()};   }
  // clang-format on

  // insert_before(pos, point) inserts the point just before pos.
  // the return value is an iterator to the newly inserted point.
  //  Note: LatticePoint && is a forwarding/universal reference, NOT a rvalue reference. Still, the
  //        function is only supposed to be called on rvalues, so we static_asserts that.
  //        This means you might have to use insert_before(pos, std::move(new_lp));


//  template <class LatticePoint, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatticePoint>>)>
//  Iterator insert_before(Iterator const &pos, LatticePoint &&new_point)
//  {
//    static_assert(std::is_reference<LatticePoint>::value == false, "Must call on rvalues");
//    return actual_list.emplace(pos.it, std::move(new_point));
//  }


  template <class LatticePoint, TEMPL_RESTRICT_DECL2(IsALatticePoint<mystd::decay_t<LatticePoint>>)>
  Iterator emplace_back(LatticePoint &&new_point)
  {
    static_assert(std::is_reference<LatticePoint>::value == false, "Must call on rvalues");
    actual_vector.emplace_back(std::move(new_point));
    auto it = actual_vector.end();
    --it;
//    return Iterator{ actual_vector.end()-- };  // this line feels soooo wrong.
    return Iterator{it};
  }

  // Inserts the lattice point pointed to by a pointer. The list takes ownership of the pointee.
  // TODO: Consider using unique_ptr
  // Untested

//  Iterator insert_before(Iterator const &pos, StoredPoint *&&point_ptr)
//  {
//    return actual_list.emplace(pos.it, std::move(point_ptr));
//  }

  // removes the element at position pos from the list and returns (and converts) it.
  // "increments" the iterator passed as argument: it now point to the next element.

  // The implementation differs, depending on whether ReturnType includes SimHashes
  // (in which case we do not recompute them)

//  template <class dummy = ReturnType, TEMPL_RESTRICT_DECL2(mystd::negation<Has_BitApprox<dummy>>)>
//  auto true_pop_point(Iterator &pos) -> ReturnType;
//
//  template <class dummy = ReturnType, TEMPL_RESTRICT_DECL2(Has_BitApprox<dummy>)>
//  auto true_pop_point(Iterator &pos) -> ReturnType;


  // erase, sort, size, empty follow std::list's semantics
  // Note that sort and size may not be available / differ for the multithreaded case.

//  Iterator erase(Iterator pos) { return actual_list.erase(pos.it); }

  void sort()
  {
    using std::sort;
//   sort(actual_vector.begin(), actual_vector.end());
  }
  typename UnderlyingContainer::size_type size() const noexcept { return actual_vector.size(); }
  NODISCARD bool empty() const noexcept { return actual_vector.empty(); }

private:
  // clang-format off
  StaticInitializer<MainStorageType> const init_stored_point_type;
//  StaticInitializer<ReturnType>  const init_return_type;
  UnderlyingContainer actual_vector;
  // clang-format on
};


// clang-format off
template<class SieveTraits>
class GaussVectorIterator<SieveTraits, false>
{
friend GaussVectorWithBitApprox<SieveTraits, false>;

private:
  using AssociatedContainerType = GaussVectorWithBitApprox<SieveTraits, false>;

//  using StoredPoint             = typename SieveTraits::GaussList_StoredPoint;
//  using ReturnType              = typename SieveTraits::GaussList_ReturnType;
  using CoordinateSelectionUsed = typename SieveTraits::CoordinateSelectionUsed;
  using MainStorageType         = typename SieveTraits::MainStorageType;

  using SimHashBlock            = typename CoordinateSelectionUsed::SimHashBlock;
  using SimHashes               = typename CoordinateSelectionUsed::SimHashes;

  // (UnderlyingIterator is only used internally (by the Iterator class in conversions)
  using UnderlyingIterator      = typename AssociatedContainerType::UnderlyingContainer::iterator;
  using CUnderlyingIterator     = typename AssociatedContainerType::UnderlyingContainer::const_iterator;
  // clang-format on

private:
  UnderlyingIterator it;

public:
  // clang-format off
  GaussVectorIterator()                                       = delete;
  GaussVectorIterator(GaussVectorIterator const &)            = default;
  GaussVectorIterator(GaussVectorIterator &&)                 = default;
  GaussVectorIterator &operator=(GaussVectorIterator const &) = default;
  GaussVectorIterator &operator=(GaussVectorIterator &&)      = default;
  // clang-format on

  // we can convert from a "plain" iterator to the underlying list. Only used internally or by
  // the list class.
private:
  // clang-format off
  constexpr GaussVectorIterator( UnderlyingIterator const &new_it) noexcept : it(new_it) {}
 // constexpr GaussIteratorBitApprox(CUnderlyingIterator const &new_it) noexcept : it(new_it) {}
  // clang-format on

public:
  // comparison, needed for for-loops (i.e. compare against list.cend() )
  bool operator==(GaussVectorIterator const &other) const { return it == other.it; }
  bool operator!=(GaussVectorIterator const &other) const { return it != other.it; }

  // increment operators. We have NO decrement, these are forward iterators.
  // (The latter is for interface compatibility with multi-threaded, where that restriction makes
  // everything both faster and easier.)
  // clang-format off
  GaussVectorIterator &operator++()    { ++it; return *this; }  // prefix version
  GaussVectorIterator  operator++(int) { return it++; }         // postfix version
  // clang-format on

  // TODO: Consider adding bool is_end() member function rather than comparing
  //       with cend() (in multithreading, cend() might mutate during an iteration,
  //       which can lead to unexpected issue.)

  // obtains an approximation to norm2 of the pointee.
  SimHashApproxNorm2 get_approx_norm2() const { return it->approx_norm2; }

  // obtains all bitapproximations of the pointee.
  SimHashes get_all_bitapproximations() const { return it->bit_approximations; }
  void update_bitapproximations() { it->update_bitapproximations(); }  // void return SimHashes?

  // gives access to the level'th block the bitapproximations.
  // In the MT case, we might return a reference to atomic.
  // NOTE: the check_simhash_scalar_product function works with *any* type that exposes an
  // access_bitapproximation(unsigned int) function (i.e. some lattice points and this class)
  // this way, we can call these directly on the iterator.
  SimHashBlock const &access_bitapproximation(unsigned int level) const
  {
    return it->bit_approximations[level];
  }

  // derefencing gives us the point stored in the STNodes (i.e. without sim_hashes)
  MainStorageType &operator*() { return *(it->ptr_to_exact); }

  // Note that for overloads of the -> operator, the return type is a class, whose -> operator is
  // recursively called (recursion stops at plain pointers).
  MainStorageType *operator->() { return it->ptr_to_exact; }

  // Iterator is explicitly convertible to pointer to const-Lattice Point.
  explicit operator MainStorageType *() const { return it->ptr_to_exact; }
};



}  // end namespace GaussSieve


#endif // include guard
