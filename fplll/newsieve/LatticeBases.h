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

namespace GaussSieve
{

/**
  LatticeBasis is the type we use internally to store lattice bases,
  along with precomputed Gram Matrix and mu matrix
  Entries:  Entries of the lattice points wrt some fixed orthonormal basis (Z^n, usually)
  gEntries: Entries of the Gram Matrix
  DimFixed: Compile-time fixed dimension (-1, if determined at runtime)
  RankFixed:Compile-time fixed rank (-1, if deteremined at runtime)

  NOTE: Entries should be a class that supports basic arithmetic operations.
  In particular, Z_NR and mpz_t will not work.

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

  constexpr PlainLatticePoint<Entries,DimFixed> get_basis_plainpoint(uint_fast16_t i) const
  {
    return make_from_any_vector<PlainLatticePoint<Entries,DimFixed>>(basis[i], ambient_dimension);
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
    if (!string_consume(is,"Dim:"))             throw bad_dumpread("LatticeBasis: Dim");
    if (!(is >> basis.ambient_dimension))       throw bad_dumpread("LatticeBasis: DimVal");
    if (!string_consume(is,"Rank:"))            throw bad_dumpread("LatticeBasis: Rank");
    if (!(is >> basis.lattice_rank))            throw bad_dumpread("LatticeBasis: RankVal");
    if (!string_consume(is,"Basis:"))           throw bad_dumpread("LatticeBasis: Basis");
    if (!read_container(is,basis.basis))        throw bad_dumpread("LatticeBasis: BasisVal");
    if (!string_consume(is,"Gram:"))            throw bad_dumpread("LatticeBasis: Gram");
    if (!read_container(is,basis.gram_matrix))  throw bad_dumpread("LatticeBasis: GramVal");
    if (!string_consume(is,"Mu:"))              throw bad_dumpread("LatticeBasis: Mu");
    if (!read_container(is,basis.mu_matrix))    throw bad_dumpread("LatticeBasis: MuVal");
    return is;
  }
};

}  // namespace GaussSieve

#endif
