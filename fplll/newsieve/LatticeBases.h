/**
  The purpose of this file is to provide glue code for the lattice basis processing capabilities
  of non-sieve fplll with the sieve part. This processing is mainly to compute GSOs.
  Some glue is needed because we do not neccessarily want to use Z_NR - types (exclusively)
  inside the sieving code, whereas the interface to compute GSOs is heavily intertwined with Z_NR.
*/

#ifndef GAUSS_SIEVE_LATTICE_BASES_H
#define GAUSS_SIEVE_LATTICE_BASES_H

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
  LatticeBasis is the type we use internally to store lattice bases,
  along with precomputed Gram Matrix and mu matrix
  Entries:  Entries of the lattice points wrt some fixed orthonormal basis (Z^n, usually)
  gEntries: Entries of the Gram Matrix
  DimFixed: Compile-time fixed dimension (-1, if determined at runtime)
  RankFixed:Compile-time fixed rank (-1, if deteremined at runtime)

  Note: We do not support a pure Gram-matrix version at the moment:
  For any Gram matrix, it is always possible to come up with some vectors that actually have
  such a gram matrix and give those to the sieve.
  Computationally, we would do that anyway, because this reduces the cost of scalar products to
  O(n) rather than O(n^2)...
  Maybe we will add such a thing for GSOLatticePoints
*/

template<class Entries, class gEntries = Entries, int DimFixed = -1, int RankFixed = -1>
class LatticeBasis
{
  static_assert(DimFixed >= -1, "Invalid DimFixed");
  static_assert(RankFixed>= -1, "Invalid RankFixed");
  static_assert( (DimFixed ==-1) || (DimFixed >= RankFixed),"Dimension must be at least rank.");
  static_assert(std::numeric_limits<uint_fast16_t>::is_specialized, "");
  static_assert(std::numeric_limits<uint_fast16_t>::max() >= DimFixed, "");

private:
  MaybeFixed<DimFixed,uint_fast16_t> ambient_dimension;

  // Technically, just number of vectors. We do not check for linear independence
  // (although the functions that construct lattice bases will probably detect this)
  MaybeFixed<RankFixed,uint_fast16_t> lattice_rank;

  // currently unused:
  // To ensure that they are >= 0
  // (in some contexts, we have the syntactic requirement are the argument is >=0, even though
  // we are guaranteed that the type / expression is never used for argument < 0 )
  // static unsigned int constexpr DimFixed_pos  = (DimFixed >= 0) ? DimFixed  : 0;
  // static unsigned int constexpr RankFixed_pos = (RankFixed>= 0) ? RankFixed : 0;

  // Triangular rank x rank Matrix (with 0 on diag)
  using MuMatrixType = ArrayOrVector< ArrayOrVector<double, RankFixed>, RankFixed>;
  MuMatrixType mu_matrix;
  // Matrix of pairwise scalar products of basis vectors
  using GramMatrixType = ArrayOrVector< ArrayOrVector<gEntries, RankFixed>, RankFixed>;
  GramMatrixType gram_matrix;
  // Basis of lattice (as a pure container object without any arithmetic defined on it)
  using BasisVector = ArrayOrVector<Entries, DimFixed>;
  using BasisType = ArrayOrVector< BasisVector, RankFixed>;
  BasisType basis;
public:
  MaybeFixed<RankFixed,uint_fast16_t> get_lattice_rank() const { return lattice_rank; }
  MaybeFixed<DimFixed,uint_fast16_t> get_ambient_dimension() const { return ambient_dimension; }
  constexpr MuMatrixType const &get_mu_matrix() const { return mu_matrix; }
  constexpr GramMatrixType const &get_g_matrix() const { return gram_matrix; }
  constexpr BasisType const &get_basis() const { return basis; }
  constexpr BasisVector const &get_basis_vector(uint_fast16_t i) const { return basis[i]; }
  constexpr Entries const &get_basis_entry(uint_fast16_t which_vec, uint_fast16_t coo) const
  {
    return basis[which_vec][coo];
  }

  double get_mu_entry(uint_fast16_t i, uint_fast16_t j) const
  {
#ifdef DEBUG_SIEVE_LOWERTRIANGULAR_MUG
    assert(j > i);
#endif
    return mu_matrix[i][j];
  }

  constexpr gEntries const &get_g_entry(uint_fast16_t i, uint_fast16_t j) const
  {
    return gram_matrix[i][j];
  }

//  constructor:
  template <class BasisContainer, class MuContainer, class GramContainer>
  LatticeBasis(uint_fast16_t const ambient_dim, uint_fast16_t const new_lattice_rank, BasisContainer &&bc, MuContainer &&muc, GramContainer &&gc)
    : ambient_dimension(ambient_dim),
      lattice_rank(new_lattice_rank),
      mu_matrix(make_array_or_vector<MuMatrixType>(muc, lattice_rank, lattice_rank)),
      gram_matrix(make_array_or_vector<GramMatrixType>(gc, lattice_rank, lattice_rank)),
      basis(make_array_or_vector<BasisType>(bc, lattice_rank, ambient_dimension))
  {
    assert(ambient_dim >= new_lattice_rank);
  }
  LatticeBasis(std::istream &is)
  {
    is >> *this;
  }

  bool operator==(LatticeBasis const &other) const
  {
    return  (ambient_dimension == other.ambient_dimension) &&
            (lattice_rank == other.lattice_rank) &&
            (mu_matrix == other.mu_matrix) &&
            (gram_matrix == other.gram_matrix) &&
            (basis == other.basis);
  }
  bool operator!=(LatticeBasis const &other) const { return !(*this == other); }

  friend std::ostream &operator<<(std::ostream &os, LatticeBasis<Entries, gEntries, DimFixed, RankFixed> const &basis)
  {
    os << "Dim: " << basis.ambient_dimension << '\n';
    os << "Rank: " << basis.lattice_rank << '\n';
    std::ios_base::fmtflags saved_flags = os.setf(std::ios_base::fixed | std::ios_base::scientific, std::ios_base::floatfield);
    os << "Basis: ";
    print_container(os,basis.basis);
    os << '\n';
    os << "Gram: ";
    print_container(os,basis.gram_matrix);
    os << '\n';
    os << "Mu: ";
    print_container(os,basis.mu_matrix);
    os << '\n';
    os.flags(saved_flags);
    return os;
  }

  friend std::istream &operator>>(std::istream &is, LatticeBasis<Entries, gEntries, DimFixed, RankFixed> &basis)
  {
    if (!string_consume(is,"Dim:")) throw bad_dumpread("LatticeBasis: Dim");
    if (!(is >> basis.ambient_dimension)) throw bad_dumpread("LatticeBasis: DimVal");
    if (!string_consume(is,"Rank:")) throw bad_dumpread("LatticeBasis: Rank");
    if (!(is >> basis.lattice_rank)) throw bad_dumpread("LatticeBasis: RankVal");
    if (!string_consume(is,"Basis:")) throw bad_dumpread("LatticeBasis: Basis");
    if (!read_container(is,basis.basis)) throw bad_dumpread("LatticeBasis: BasisVal");
    if (!string_consume(is,"Gram:")) throw bad_dumpread("LatticeBasis: Gram");
    if (!read_container(is,basis.gram_matrix)) throw bad_dumpread("LatticeBasis: GramVal");
    if (!string_consume(is,"Mu:")) throw bad_dumpread("LatticeBasis: Mu");
    if (!read_container(is,basis.mu_matrix)) bad_dumpread("LatticeBasis:MuVal");
    return is;
//  if basis.read_dim(is) == false throw bad_dumread
  }
};

/**
  SieveLatticeBasis is the type used to represent lattice bases.
  A basis of this type is stored inside the sieve.
  We support converting to an array of PlainLatticePoints and extracting GSO information.

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
// clang-format off
template<class SieveTraits, bool MT, class Enabled>
class SieveLatticeBasis
// clang-format on
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

  // Note : InputET_NOZNR is not Z_NR - wrapped (because of the way ZZMat works...)
  using InputET_NOZNR      = typename IsZZMatClass<InputBasisType>::GetET;
  // InputET is the Z_NR - type compatible with the ZZMat-type of the input basis.
  using InputET            = fplll::Z_NR<InputET_NOZNR>;

  // Same as ET_NOZNR, but ET_NOZNRFixed may be mpz_class instead of mpz_t
  using InputET_NOZNRFixed = PossiblyMpztToMpzClass<InputET_NOZNR>;

  // TODO: This makes little sense... LengthType has a different meaning.
  using OutputET           = typename SieveTraits::LengthType;
  // clang-format on

  using DimensionType               = typename SieveTraits::DimensionType;
  using GlobalStaticDataInitializer = typename SieveTraits::GlobalStaticDataInitializer;

  // This one is hard-coded to a PlainLatticePoint
  using BasisVectorType = PlainLatticePoint<OutputET, SieveTraits::get_nfixed>;

  using GSOType = fplll::MatGSO<InputET, fplll::FP_NR<double>>;

public:
  // Note: We copy the basis into originial_basis.
  // This is because the GSO object actually uses a reference and modifies original_basis.
  // clang-format off
  explicit SieveLatticeBasis(InputBasisType const &input_basis,
                             GlobalStaticDataInitializer const &static_data)
      : original_basis(input_basis),
        ambient_dimension(input_basis.get_cols()),
        init_basis_vector_type(static_data),
        lattice_rank(input_basis.get_rows()),
        // u,u_inv intentionally uninitialized
        // will be initialized in call to
        // GSO(original_basis, u, u_inv, fplll::MatGSOInterfaceFlags::GSO_INT_GRAM),

        // mu_matrix and g_matrix get initialized as empty rank x rank matrices
        mu_matrix(lattice_rank, std::vector<double>(lattice_rank)),
        g_matrix(lattice_rank, std::vector<InputET_NOZNRFixed>(lattice_rank)),
        basis_vectors()

  // clang-format on
  {
    fplll::Matrix<InputET> u, u_inv;  //, g;
    GSOType GSO(original_basis, u, u_inv, fplll::MatGSOInterfaceFlags::GSO_INT_GRAM);
    GSO.update_gso();  // todo: raise exception in case of error
    // We compute these at creation to simplify thread-safety issues.
    compute_mu_matrix(GSO);
    compute_g_matrix(GSO);

    // extract and convert the actual lattice vectors.
    basis_vectors.reserve(lattice_rank);
//    basis_vectors = new BasisVectorType[lattice_rank];
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      // make_from_znr_vector copies and push_back uses move semantics.
      basis_vectors.push_back( make_from_znr_vector<BasisVectorType>(input_basis[i], ambient_dimension) );
    }
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
  void compute_mu_matrix(GSOType &GSO)  // Note: GSO is not passed as const-ref.
  {
    // use GSO's capabilities and convert to non-FP_NR - type
    fplll::Matrix<fplll::FP_NR<double>> ZNR_mu = GSO.get_mu_matrix();
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      for (uint_fast16_t j = 0; j < lattice_rank; ++j)
      {
        mu_matrix[i][j] = ZNR_mu[i][j].get_d();
      }
    }
    return;
  }

  void compute_g_matrix(GSOType &GSO)  // Note: GSO is not passed as const-ref.
  {
    fplll::Matrix<InputET> const gmatrix_GSO = GSO.get_g_matrix();  // class returned by GSO
    // convert fplll::Matrix<ET> to vector<vector<ET_NOZNRFixed>>
    for (uint_fast16_t i = 0; i < lattice_rank; ++i)
    {
      for (uint_fast16_t j = 0; j < lattice_rank; ++j)
      {
        g_matrix[i][j] = static_cast<InputET_NOZNRFixed>(gmatrix_GSO(i, j).get_data());
//        std::cout << g_matrix[i][j] << " ";
      }
//      std::cout << std::endl;
    }
  }

  std::vector<std::vector<double>> get_mu_matrix() const { return mu_matrix; }

  double get_maxbistar2() const { return maxbistar2; }

  std::vector<std::vector<InputET_NOZNRFixed>> get_g_matrix() const { return g_matrix; }

  // returns (i,j)th entry of mu. Assumes j>i
  double get_mu(uint_fast16_t i, uint_fast16_t j) const
  {
#ifdef DEBUG_SIEVE_LOWERTRIANGULAR_MUG
    assert(j > i);
#endif
    return (mu_matrix)[i][j];
  }

  // returns (i,j)th entry of g. Assumes j>=i
  InputET_NOZNRFixed get_g(uint_fast16_t i, uint_fast16_t j) const
  {
#ifdef DEBUG_SIEVE_LOWERTRIANGULAR_MUG
    assert(j >= i);
#endif  // DEBUG_SIEVE_LOWERTRIANGULAR_MUG
    return (g_matrix)[i][j];
  }

  // Returns ith basis vector.
  BasisVectorType const &get_basis_vector(uint_fast16_t i) const
  {
    assert(i < lattice_rank);
    return basis_vectors[i];
  }
  // we have  1/(2*Pi*exp(1)) =  0.05854983150 for GH
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
    mink_bound          = static_cast<InputET_NOZNRFixed>(mink_bound_d);
    //std::cout << "root_det = " << root_det  << std::endl;
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
#endif

  InputET_NOZNRFixed get_minkowski_bound() const { return mink_bound; }

private:
  InputBasisType original_basis;

public:
  DimensionType const ambient_dimension;
  StaticInitializer<BasisVectorType> const init_basis_vector_type;
  uint_fast16_t const lattice_rank;  // Technically, just number of vectors.
                                     // We don't verify linear independence ourselves.
                                     // (even though GSO computation does, probably)
#ifdef PROGRESSIVE
  std::vector<double> progressive_bounds; // progressive_bounds[i] stores GH^2 for
                                          // L[n_start+i] (to determine when we increase the rank)
#endif
private:
  //  fplll::Matrix<InputET> u, u_inv; //, g;
  //  fplll::MatGSO<InputET, fplll::FP_NR<double>> GSO;
  // precomputed on demand:
  // clang-format off
  std::vector<std::vector<double            >> mu_matrix;
  std::vector<std::vector<InputET_NOZNRFixed>> g_matrix;
  // clang-format on
  InputET_NOZNRFixed mink_bound;

  std::vector<BasisVectorType> basis_vectors;
  double maxbistar2;
};  // end of class

}  // namespace GaussSieve

#endif
