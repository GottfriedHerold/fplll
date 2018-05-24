#ifndef LATTICE_BASIS_CONVERT_H
#define LATTICE_BASIS_CONVERT_H

#include "DefaultIncludes.h"
#include "PlainLatticePoint.h"
#include "SieveUtility.h"
#include "Typedefs.h"
#include "fplll/defs.h"
#include "fplll/gso.h"
#include "fplll/nr/matrix.h"
#include "fplll/nr/nr_Z.inl"
#include "gmp.h"
#include "gmpxx.h"

namespace GaussSieve
{

/**
  SieveLatticeBasis is used to convert whatever input we have into a LatticeBasis.
  The purpose of this class is to have all the Z_NR - ugliness in one place for the sieving code.
  SieveLatticeBasis holds a (smart) pointer to a LatticeBasis

  Notably SieveLatticeBasis(InputBasisType const &input_basis, [...] ) takes
  an input_basis (possibly of Z_NR - type) and converts it into the internal types that we use.
  The allowed InputBasisTypes supported are controlled by SieveTraits.
*/

// NOTE: The last template parameter Enabled is a dummy parameter and is always void.
// Its purpose is to allow partial specializations that depend on properties of SieveTraits.
//
// Notably, you can do partial specializations
// template<class SieveTraits, bool MT>
// class SieveLatticeBasis<SieveTraits, MT, mystd::enable_if_t< CONDITION > >
// {
// ...
// };
// where CONDITION is a bool that depends on SieveTraits and MT.
// Since template argument deduction takes places before template argument substitution and SFINAE,
// this has the effect that the specialization is only considered iff CONDITION is true.

// This way, we can write various versions, depending on whether we receive Z_NR - input types or
// not. Currently, we only take ZZ_Mat - classes

template <class SieveTraits, bool MT, class Enabled = void> class SieveLatticeBasis;

// unusable default (by design)
template<class SieveTraits, bool MT, class Enabled>
class SieveLatticeBasis
{
  // This must never be instantiated.
  static_assert(std::is_void<Enabled>::value == false, "");  // this will always trigger!
  static_assert(SieveTraits::IsSieveTraitsClass::value, "Invalid Traits class");
public:
  SieveLatticeBasis(...) = delete;
};

/**
  Specialization for classes where SieveTraits::InputBasisType is recognized as a ZZMat-class.
  This is the only implementation we currently have.
  This essentially wraps around the fplll GSO class interfac(es).

  Note that the GSO class interfaces may not have been designed with multi-threading in mind.
  In particular, they perfom lazy evaluation and store values already computed.
  So even what should be read-only operations are probably not thread-safe.
  For that reason, we only use the GSO methods during creation of the SieveLatticeBasis object and
  store the results in our own classes. This is quite suboptimal and ugly, but it's not dominating
  the cost of the whole algorithm anyway.

  TODO: Specify interface
*/

// TODO: Make GSO object local to constructor.

// clang-format off
template <class SieveTraits, bool MT>
class SieveLatticeBasis<SieveTraits, MT,  // next param selectively enables this specialization:
    mystd::enable_if_t<IsZZMatClass<typename SieveTraits::InputBasisType>::value> >
// clang-format on
{

private:
  // clang-format off
  // Transform the received InputBasisType into something appropriate...
  using InputBasisType     = typename SieveTraits::InputBasisType;
  static_assert(IsZZMatClass<InputBasisType>::value,""); // checked by enable_if anyway.

  // Get template arg: InputBasisType = ZZMat<InputET_NOZNR>
  // Note : InputET_NOZNR is not Z_NR - wrapped (because of the way ZZMat works...), hence the name
  using InputET_NOZNR      = typename IsZZMatClass<InputBasisType>::GetET;
  // InputET is the Z_NR - type compatible with the ZZMat-type of the input basis.
  using InputET_ZNR        = fplll::Z_NR<InputET_NOZNR>;

  // Same as ET_NOZNR, but ET_NOZNRFixed may be mpz_class instead of mpz_t
  using InputET_Fixed = PossiblyMpztToMpzClass<InputET_NOZNR>;

  // TODO: This makes little sense... LengthType has a different meaning.
//  using OutputET           = InputET_NOZNRFixed;
  using OutputET           = typename SieveTraits::BasisEntries;
  using LengthType         = typename SieveTraits::LengthType;
  // clang-format on

  using DimensionType               = typename SieveTraits::DimensionType;
  using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;

  // This one is hard-coded to a PlainLatticePoint
//  using BasisVectorType = PlainLatticePoint<OutputET, SieveTraits::get_nfixed>;
  using BasisVectorType = typename SieveTraits::PlainPoint;
  using OutputBasis = typename SieveTraits::BasisType;

  using GSOType = fplll::MatGSO<InputET_ZNR, fplll::FP_NR<double>>;

private:
//  InputBasisType original_basis;

public:
  DimensionType const ambient_dimension;
  StaticInitializer<BasisVectorType> const init_basis_vector_type;
  uint_fast16_t const lattice_rank;  // Technically, just number of vectors.
                                     // We don't verify linear independence ourselves.
                                     // (even though GSO computation does, probably)
private:
  //  fplll::Matrix<InputET> u, u_inv; //, g;
  //  fplll::MatGSO<InputET, fplll::FP_NR<double>> GSO;
  // precomputed on demand:
  // clang-format off
//  std::vector<std::vector<double            >> mu_matrix;
//  std::vector<std::vector<InputET_NOZNRFixed>> g_matrix;
  // clang-format on
  LengthType mink_bound;

  std::vector<BasisVectorType> copy_basis_vectors;
  double maxbistar2;

//  using OutputBasis = LatticeBasis<OutputET, LengthType, SieveTraits::get_nfixed, SieveTraits::get_nfixed>;

  std::shared_ptr<OutputBasis> ptr_basis;
#ifdef PROGRESSIVE
  std::vector<double> progressive_bounds; // progressive_bounds[i] stores GH^2 for
                                           // L[n_start+i] (to determine when we increase the rank)
#endif


public:

  OutputBasis const &get_basis() const { return *ptr_basis; }
  std::shared_ptr<OutputBasis> get_basis_ptr() const { return ptr_basis; }

  // Note: We copy the basis into originial_basis.
  // This is because the GSO object actually uses a reference and modifies original_basis.
  // clang-format off
  explicit SieveLatticeBasis(InputBasisType const &input_basis,
                             GlobalStaticDataInitializer const &static_data)
      :
//       original_basis(input_basis),
        ambient_dimension(input_basis.get_cols()),
        init_basis_vector_type(static_data),
        lattice_rank(input_basis.get_rows()),
        // mu_matrix and g_matrix get initialized as empty rank x rank matrices
//        mu_matrix(lattice_rank, std::vector<double>(lattice_rank)),
//        g_matrix(lattice_rank, std::vector<InputET_NOZNRFixed>(lattice_rank)),
        copy_basis_vectors()

  // clang-format on
  {
    InputBasisType original_basis ( input_basis );  // copy input basis, because GSO object
                                                    // modifies it.
    fplll::Matrix<InputET_ZNR> u, u_inv;
    GSOType GSO(original_basis, u, u_inv, fplll::MatGSOInterfaceFlags::GSO_INT_GRAM);
    GSO.update_gso();  // todo: raise exception in case of error
    // create empty rank x rank matrices:
    std::vector<std::vector<double>> mu_matrix ( lattice_rank, std::vector<double>(lattice_rank));
    std::vector<std::vector<double>> r_matrix ( lattice_rank, std::vector<double>(lattice_rank));
    std::vector<std::vector<LengthType>> g_matrix ( lattice_rank, std::vector<LengthType>(lattice_rank));
    std::vector<std::vector<OutputET>> basis_entries ( lattice_rank, std::vector<OutputET>(ambient_dimension));
    // We compute these at creation to simplify thread-safety issues.
    std::cout << "XXX" << '\n';
    compute_mu_matrix(GSO, mu_matrix, lattice_rank);
    print_container(std::cout, mu_matrix);

    std::cout << "XXX" << '\n';
    compute_r_matrix(GSO, r_matrix, lattice_rank);
    print_container(std::cout, r_matrix);

    std::cout << "XXX" << '\n';
    compute_g_matrix(GSO, g_matrix, lattice_rank);
    print_container(std::cout, g_matrix);
    std::cout << "XXX" << '\n';

    // extract and convert the actual lattice vectors.
//    basis_vectors.reserve(lattice_rank);
//    basis_vectors = new BasisVectorType[lattice_rank];
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      for (uint_fast16_t j = 0; j < ambient_dimension; ++j)
      {
        basis_entries[i][j] = convert_to_inttype<OutputET>(input_basis[i][j].get_data());
      }
//      convert_to_inttype
//      // make_from_znr_vector copies and push_back uses move semantics.
      copy_basis_vectors.push_back( make_from_znr_vector<BasisVectorType>(input_basis[i], ambient_dimension) );
    }
    ptr_basis = std::make_shared<OutputBasis>(ambient_dimension, lattice_rank, basis_entries, mu_matrix, r_matrix, g_matrix);
    maxbistar2 = GSO.get_max_bstar().get_d();

    compute_minkowski_bound(GSO);

#ifdef PROGRESSIVE
    progressive_bounds.resize(lattice_rank);
    compute_progressive_bounds(GSO);
#endif
  }

  // clang-format off
  SieveLatticeBasis(SieveLatticeBasis const &)            = delete;
  SieveLatticeBasis(SieveLatticeBasis &&)                 = default;
  SieveLatticeBasis &operator=(SieveLatticeBasis const &) = delete;
  SieveLatticeBasis &operator=(SieveLatticeBasis &&)      = default;
  // clang-format on

  ~SieveLatticeBasis() = default; // { delete[] basis_vectors; }

  // Note: Const-correctness is strange wrt. the fplll::GSO classes.
  // Consider "patching" this up by marking GSO mutable.
  // TODO: Patch the MatGSO class upstream instead (use mutable for lazy evaluation...)

  // IMPORTANT: get_mu_matrix()[i][j] is only meaningful for j>i!

  // helper function that precomputes, converts and stores the whole mu_matrix:
private:
  // Note: GSO is not passed as const-ref.
  template<class Container>
  static void compute_mu_matrix(GSOType &GSO, Container &mu_matrix, uint_fast16_t const lattice_rank)
  {
    // use GSO's capabilities and convert to non-FP_NR - type
    fplll::Matrix<fplll::FP_NR<double>> ZNR_mu = GSO.get_mu_matrix();
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      for (uint_fast16_t j = 0; j < lattice_rank; ++j)
      {
        if (j < i)
        {
          mu_matrix[i][j] = ZNR_mu[i][j].get_d();
        }
        else
        {
          mu_matrix[i][j] = 0;
        }
      }
    }
    return;
  }

  template<class Container>
  static void compute_r_matrix(GSOType &GSO, Container &r_matrix, uint_fast16_t const lattice_rank)
  {
    fplll::Matrix<fplll::FP_NR<double>> ZNR_r = GSO.get_r_matrix();
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      for (uint_fast16_t j = 0; j < lattice_rank; ++j)
      {
        if (j < i)
        {
          r_matrix[i][j] = ZNR_r[i][j].get_d();
        }
        else
        {
          r_matrix[i][j] = 0;
        }
      }
    }
    return;
  }

  // Note: GSO is not passed as const-ref.
  template<class Container>
  static void compute_g_matrix(GSOType &GSO, Container &g_matrix, uint_fast16_t const lattice_rank)
  {
    fplll::Matrix<InputET_ZNR> const gmatrix_GSO = GSO.get_g_matrix();  // class returned by GSO
    // convert fplll::Matrix<ET> to vector<vector<ET_NOZNRFixed>>
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      for (uint_fast16_t j = 0; j < lattice_rank; ++j)
      {
        g_matrix[i][j] = convert_to_inttype<LengthType>(gmatrix_GSO(i, j).get_data());
      }
    }
    // make symmetric: fplll's GSO get_g_matrix might (depending on version) only compute and access
    // entries with i <= j.
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      for (uint_fast16_t j = 0; j < lattice_rank; ++j)
      {
        if(i < j)
        {
          g_matrix[i][j] = g_matrix[j][i];
        }
      }
    }
  }
public:

//  std::vector<std::vector<double>> get_mu_matrix() const { return ptr_basis->get_mu_matrix; }

  double get_maxbistar2() const { return maxbistar2; }

//  std::vector<std::vector<InputET_NOZNRFixed>> get_g_matrix() const { return ptr_basis -> get_g_matrix; }

  // returns (i,j)th entry of mu. Assumes j>i
  double get_mu(uint_fast16_t i, uint_fast16_t j) const
  {
    return ptr_basis->get_mu_entry(i,j);
  }

  // returns (i,j)th entry of g. Assumes j>=i
  InputET_Fixed get_g(uint_fast16_t i, uint_fast16_t j) const
  {
    return ptr_basis->get_g_entry(i,j);
  }

  // Returns ith basis vector.
  BasisVectorType const &get_basis_vector(uint_fast16_t i) const
  {
    assert(i < lattice_rank);
    return copy_basis_vectors[i];
  }
//  // we have  1/(2*Pi*exp(1)) =  0.05854983150 for GH
  // 0.644 is ok

  void compute_minkowski_bound(GSOType &GSO)
  {
    static constexpr double gh_prefactor = 0.592399062725144101295483424139;
    // 1.05 / sqrt(pi)
    static constexpr double gh_prefactor2 = gh_prefactor * gh_prefactor;

    // returns det(B)^{2/dim}
    // fplll::FP_NR<double> root_det2 =  GSO.get_root_det(1, lattice_rank);
    double root_det     = GSO.get_root_det(0, lattice_rank).get_d();
    double mink_bound_d = gh_prefactor2 * std::pow( std::tgamma(0.5 * lattice_rank + 1), 2. / lattice_rank)
                          * root_det;
//    double mink_bound_d = 0.0644*root_det * static_cast<double>(lattice_rank);
    mink_bound          = static_cast<LengthType>(mink_bound_d);
    std::cout << "mink_bound is set to: " << mink_bound << std::endl;
  }

#ifdef PROGRESSIVE
  void compute_progressive_bounds(GSOType &GSO)
  {
    for (unsigned int i = 0; i<lattice_rank; ++i)
    {
      // when progressive_rank == i+1;
      // GSO.get_root_det (supposedly) returns the determinant of (i+1)-dim. sublattice
      progressive_bounds[i] =0.0664 * (i+1) * convert_to_double(  GSO.get_root_det(0, i+1).get_d() );
      //std::cout << "progressive_bounds[i] = " << progressive_bounds[i] << std::endl;
    }
  }

  double get_progressive_bound (uint_fast16_t const i) const
  {
    return progressive_bounds[i];
  }

#endif

  LengthType get_minkowski_bound() const { return mink_bound; }
};  // end of class


} // end namespace GaussSieve

#endif // include guard
