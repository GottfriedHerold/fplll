#pragma once
#ifndef GAUSS_SIEVE_BITAPPROX_NEW_H
#error wrong usage
#endif

namespace GaussSieve
{

/**
  Constructor of BlockOrthogonalSimHash from a dimension and a random seed
*/

template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
BlockOrthogonalSimHash<sim_hash_len,sim_hash_num,MT,DimensionType_arg>::
    BlockOrthogonalSimHash(DimensionType const dim, unsigned int random_seed)
{
  std::mt19937 rng;  // Mersenne twister random engine. We don't really care about quality anyway.
  rng.seed(random_seed);  // seed with the input random seed

  unsigned int const ambient_dimension = dim;  // Note that this converts dim to an unsigned int.
  assert(ambient_dimension!=0);

  // sim_hash_num * sim_hash_len is the total number of required bits in the output simhashes.
  // We need ceil(total number of required bits / ambient_dimension) many repetitions of the input
  // vector to get at least that many output bits.
  // Observe that  ((X-1) / Y) + 1 is the same as ceil(X/Y) for unsigned integers X,Y
  number_of_orthogonal_blocks = ((sim_hash_num * sim_hash_len - 1) / ambient_dimension) + 1;

  // fast_walsh_hadamard_len is the largest power of 2 which is <= ambient dimension.
  fast_walsh_hadamard_loglen  = static_cast<unsigned int>(floor(log2(ambient_dimension)));
  fast_walsh_hadamard_len     = static_cast<unsigned int>(pow(2, floor(log2(ambient_dimension))));

  // pmatrices are a std::vector<std::array<pmatrix<num_of_transforms>>.
  // The outer dimension (of the vector) of given by number_of_orthogonal_blocks. Dito dmatrices.
  pmatrices.resize(number_of_orthogonal_blocks);
  dmatrices.resize(number_of_orthogonal_blocks);
  // construct each pmatrix / dmatrix using the corresponding constructor.
  // This means all patrices / dmatrices are independent. (and uniform, at the moment)
  for (unsigned int i = 0; i < number_of_orthogonal_blocks; ++i)
  {
    for (unsigned int j = 0; j < num_of_transforms; ++j)
    {
      pmatrices[i][j] = PMatrix(ambient_dimension,rng);
      dmatrices[i][j] = DMatrix(ambient_dimension,rng);
    }
  }
}


/**
 fast_partial_walsh_hadamard<len>(input) performs (fast) Walsh-Hadamard Transform on the first
 len coordinates of its input. len is enforced to be a power of two.
 The entries of the input vector must support arithmentic (+,-) and be default-constructible and
 swappable.

 NOTE: This function is deprecated and not maintained or used.
*/

// Note: input is passed by-value. We modify the local copy.
template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
template<class T>
inline auto BlockOrthogonalSimHash<sim_hash_len,sim_hash_num,MT,DimensionType_arg>::
    fast_partial_walsh_hadamard(std::vector<T> input) const
    -> std::vector<T>
{
  unsigned int const len = fast_walsh_hadamard_len;
  using std::swap;
  assert(is_a_power_of_two(len));  // len must be a power of two:
  assert(len <= input.size());
  std::vector<T> output(input.size());  // will eventually hold the output

  for (uint_fast16_t i = len; i < output.size(); ++i)  // copy coordinates untouched by the trafo.
  {
    output[i] = input[i];
  }

  // on [0...len-1] coordinates perform WH

  for (uint_fast16_t i = len >> 1; i > 0; i >>= 1)
  {
    for (uint_fast16_t j = 0; j < len; j++)
    {
      // Note: static_cast<T> is required to be able to use mpz_class, because mpz_class internally
      // uses expression templates (so the types of a+b and a-b don't match for mpz_classes a,b)
      output[j] = ( (j/i)%2 != 0) ? static_cast<T>(input[j-i] - input[j]  )   // for j/i odd
                                  : static_cast<T>(input[j]   + input[i+j]);  // for j/i even
    }
    swap(input, output);
  }

  const double lengthfactor = std::sqrt(len);  // we have to properly rescale the modified coos.
  for(uint_fast16_t i = 0; i < len; ++i)
  {
    output[i] /= lengthfactor;
  }
  return output;
}

/**
  faster_almost_partial_walsh_hadamard_inplce(input) performs an inplace (i.e. changing input)
  (almost) WH-transform on the first fast_walsh_hadamard_len coordinates.
  (fast_walsh_hadamard_len is a member field of the class and should be the largest power of 2 that
   fits into the dimension of the vectors on which we operate)
  The "almost" refers to the fact that we actually perform the WH-trafo up to a fixed permutation of
  coos and sign changes (permutations will not affect which coos are affected, just permute among
  those affected).
  This version is much faster than the above one.
*/

template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
template<class T>
void inline BlockOrthogonalSimHash<sim_hash_len, sim_hash_num, MT, DimensionType_arg>::
    faster_almost_partial_walsh_hadamard_inplace(std::vector<T> & input) const
{
  unsigned int const len = fast_walsh_hadamard_len;
  assert(is_a_power_of_two(len));
  assert(len <= input.size());

  T tmp;

  // Note: The (fast) WH transform can be viewed as a sequence of log2(len) matrices A_i, where each
  // matrix is (after a permutation of coordinates) a block diagonal matrix, where each block is
  // the 2x2 matrix (up to scaling)
  //           (1   1)
  // A_basic = (     )  (up to scaling by sqrt(2) )
  //           (1  -1)

  // The permutations depend on i; notably for A_i, the blocks correspond to the pair of indices
  // that differ only in the i'th bit (counting from msb)
  // To perform a fast WH trafo, we just apply a (variant) of A_basic to such a pair of coos.
  // Note that the key difference to the version above is that applying such a 2x2 matrix only
  // uses ONE temporary (scalar) variable, as opposed to a temporary matrix. The lower memory foot-
  // print is what gains an order of magnitude in speed.



  // writing the indices in binary, for the WH-Transform, we perform
  // output[****0*****] = input[*****0*****] + input[*****1*****]
  // output[****1*****] = input[*****0*****] - input[*****1*****]
  // and then set input to output, where the affected "special" bit ranges over all possibilites,
  // starting with the msb.

  // Changing the processing order for the "special bit" and where the -sign appears above does
  // not matter for our application:
  // It is equivalent to a fixed permutation of coos and multiplying some coos by +/-1 before/after
  // the WH transformation; choose a variant that requires slightly less computation, notably we can
  // use output[****1****] = -input[*****0*****] + input[****1****] to use a single -= operation
  // (this is faster by ~5 -- 10%)


  // The following 3 versions are equivalent, up to permutation and signs
  // experimentally, the one selected was the fastest on my (Gotti's) machine...

  // Version 1:

  // i            = 0000010000 in the above notation
  // higher_bits is *****00000
  // lower_bits  is 000000****
  // (i.e. the loop variables are already multiplied by the appropriate power of 2)

  for (uint_fast16_t i = len >> 1; i > 0; i >>= 1)
//  for (uint_fast16_t i = 1; i < len; i<<=1) // loop order does not matter (up to sign/permutation)
  {
    // this order of the 2 inner loops is 5-10% faster. No idea why.
    for(uint_fast16_t lower_bits = 0; lower_bits < i; ++lower_bits)
    {
      for(uint_fast16_t higher_bits =0; higher_bits < len; higher_bits += (2*i) )
      {
        // Note all + in the indices are the same as XOR's
        tmp = input[higher_bits+0+lower_bits];
        input[higher_bits+0+lower_bits]+=input[higher_bits+i+lower_bits];
//       input[higher_bits+i+lower_bits] =tmp - input[higher_bits+i+lower_bits];   // <- WH trafo
        input[higher_bits+i+lower_bits]-=tmp;  //  <- 5-10% faster.
      }
    }
  }

  /*
  // Version 2 : Do not pre-multiply the loop variables by the corresponding powers of 2 and instead
  //             bit-shift them when needed.
  unsigned int const log_len = fast_walsh_hadamard_loglen;
  for(uint_fast16_t i=0; i < log_len;++i)
  {
    unsigned int const ip = (1 << i);

    for(uint_fast16_t j=0; j < (1 << i); ++j)
    {
      for(uint_fast16_t k=0; k < (1 << (log_len-1-i)); ++k)
      {
        tmp = input[j+ (k<<(i+1))];
        input[j+(k<<(i+1))] += input[j+(k<<(i+1))+ip];
        input[j+(k<<(i+1))+ip]-=tmp;
      }
    }
  }
  */

  /*
  // Version 3:
  // We iterate over all indices j that have a 0 on the bitposition marked by i by testing the bit.
  for (uint_fast16_t i = len >> 1; i > 0; i >>= 1)
  {
    for(uint_fast16_t j = 0; j < len; ++j)
    {
      if ((j & i) !=0) continue;
      tmp = input[j];
      input[j]+=input[j^i];
      input[j^i]-=tmp;
    }
  }
  */

  double const lengthfactor = std::sqrt(len);  // we have to properly rescale the modified coos.
  for(uint_fast16_t i = 0; i < len; ++i)
  {
    input[i] /= lengthfactor;
  }
  return;
}

/**
 computes sim_hash of a lattice point
*/

template<std::size_t sim_hash_len, std::size_t sim_hash_num, bool MT, class DimensionType_arg>
template<class LatP, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatP>)>
inline auto BlockOrthogonalSimHash<sim_hash_len,sim_hash_num,MT,DimensionType_arg>::
    compute_all_bitapproximations(LatP const &point) const -> SimHashes
{
  unsigned int const dim = static_cast<unsigned int>(point.get_dim());
  // Note: we currently do not check whether this dim is the same as the dim used to initialize
  //       the class (which we assume to be the case)
  // TODO: Check that, at least in debug mode
  using ET    = Get_AbsoluteCooType<LatP>;
  using Block = std::vector<ET>;  // will be of length dim.
  // TODO: Optimize for fixed dimension.

  // copy the input point into a std::vector<ET>.
  // This is needed, because we no not know the internal details of the input lattice point.
  // in particular, we neither know that the lattice point internally directly stores the data we
  // need do not want to maintain any extra data from the lattice point (such as its length)
  // TODO: Maybe use PlainLatticePoint?
  Block copy_of_point(dim);
  for (unsigned int i = 0; i < dim; ++i)
  {
    copy_of_point[i] = point.get_absolute_coo(i);
  }

  // initialize blocks by number_of_orthogonal_blocks many copies of the input point.
  std::vector<Block> blocks(number_of_orthogonal_blocks, std::move(copy_of_point));

  // each such block gets transformed by applying WH * D * P num_transform many times.
  for (uint_fast16_t i = 0; i < number_of_orthogonal_blocks; ++i)  // for each block:
  {
    for (uint_fast8_t j = 0; j < num_of_transforms; ++j)  // repeat num_of_transforms many times:
    {
      pmatrices[i][j].apply(blocks[i]);
      dmatrices[i][j].apply(blocks[i]);
      if ( j < num_of_transforms - 1 )  // on the last iteration, we do not perform a WH - trafo
      {
        faster_almost_partial_walsh_hadamard_inplace(blocks[i]);
      }
    }
  }

  // put together the blocks into an array of bitsets.
  SimHashes ret;
  for (unsigned int n = 0; n < sim_hash_num; ++n)
  {
    for (unsigned int m = 0; m < sim_hash_len; ++m)
    {
      // index of the bit currently considered
      // if we use only one level of indexing.
      unsigned int const flat_bit_count = n * sim_hash_len + m;

      ret[n][m] = (blocks[flat_bit_count / dim][flat_bit_count % dim] > 0);
    }
  }
  return ret;
}

/**
  construct a random permutation in dim dimension. Just uses std::shuffle for that.
*/
PMatrix::PMatrix(unsigned int dim, std::mt19937 &rng)
{
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("about to initialize P matrix")
  permutation.resize(dim);
  for (uint_fast16_t i = 0; i < dim; ++i)
  {
    permutation[i] = i;
  }
  // TODO:  Consider making the permutation fixed / improve the mixing between the parts affected by
  //        WH-Trafo. With the current amount of randomness, some coordinates are unaffected.
  std::shuffle(permutation.begin(), permutation.end(), rng);
#ifdef DEBUG_SIEVE_TRACEINITIATLIZATIONS
//  print();
#endif
}

/**
  Applies stored permutation to a vector.
*/
template<class T>
inline void PMatrix::apply(std::vector<T> &vec) const
{
  auto const dim = vec.size();
  assert(dim == permutation.size());
  std::vector<T> out(dim);  // temporary holding the new vec after the function finishes
  for (unsigned int i = 0; i < dim; ++i)
  {
    out[i] = vec[permutation[i]];
  }
  using std::swap;
  swap(out, vec);
  // Note: We could avoid the temporary by storing a cycle representation of the permutation.
  //       At any rate, we might remove the randomness from the permutation...
  return;
}

/**
  Print the stored data to stream. Only used for debugging at the moment.
  TODO: Turn into << operator
*/
inline void PMatrix::print(std::ostream &os) const
{
  // os << "P  [" << i <<"] is: " << std::endl;
  for (uint_fast16_t k = 0; k < permutation.size(); ++k)
  {
    os << permutation[k] << " ";
  }
  os << std::endl;
}

/**
  Construct a (uniformly random) diagonal matrix with +/-1 on the diagonals.
*/

DMatrix::DMatrix(unsigned int const dim, std::mt19937 &rng)
{
  DEBUG_SIEVE_TRACEINITIATLIZATIONS("about to fill-up the D matrix")
  std::uniform_int_distribution<uint_fast8_t> distr(0, 1);
  diagonal.resize(dim);
  for (uint_fast16_t i = 0; i < dim; ++i)
  {
    diagonal[i] = distr(rng);
  }
#ifdef DEBUG_SIEVE_TRACEINITIATLIZATIONS
//  print();
#endif
}

/**
  Applies such a diagonal matrix on a vector (i.e. flip some signs in vec).
*/
template<class T> inline void DMatrix::apply(std::vector<T> &vec) const
{
  auto const dim = vec.size();
  assert(dim == diagonal.size());
  for (unsigned int i = 0; i < dim; ++i)
  {
    if (diagonal[i] != 0)  // otherwise, keep vec[i] as it is.
    {
      vec[i] = -vec[i];
    }
  }
}

/**
  Print the stored data to stream. Only used for debugging at the moment.
  TODO: Turn into << operator
*/
inline void DMatrix::print(std::ostream &os) const
{
  for (uint_fast16_t k = 0; k < diagonal.size(); ++k)
  {
    os << static_cast<int>(diagonal[k]) << " ";
  }
  os << std::endl;
}

}  // end namespace GaussSieve


/**
  UNUSED CODE STARTS HERE:
*/


/*
// unused, did not incorporate recent changes.
// (unused) variant of the above function, for arrays instead of vectors.
template<class T, std::size_t arraylen>
inline auto fast_partial_walsh_hadamard(std::array<T,arraylen> input, unsigned int const len)
    -> std::array<T,arraylen>
{
  using std::swap;
  assert(std::bitset< std::numeric_limits<unsigned long>::digits>{len}.count() == 1);
  assert(len <= arraylen);
  std::array<T,arraylen> output{}; // assumes that entries are default-constructible.

  const double lengthfactor = std::sqrt(len);  // we have to properly rescale the unmodified coos.
  for (uint_fast16_t i = len; i < arraylen; ++i)
  {
    input [i] = static_cast<T>(input[i] * lengthfactor);
    output[i] = input[i];
  }

  for (uint_fast16_t i = len >> 1; i > 0; i >>= 1)
  {
    for (uint_fast16_t j = 0; j < len; j++)
    {
      // Note: static_cast<T> is required to be able to use mpz_class, because mpz_class internally
      // uses expression templates (so the types of a+b and a-b don't match for mpz_classes a,b)
      output[j] = ( (j/i)%2 != 0) ? static_cast<T>(input[j-i] - input[j]  )   // for j/i odd
                                  : static_cast<T>(input[j]   + input[i+j]);  // for j/i even
    }
    swap(input, output);
  }
  return output;
}
*/

#if 0  // deactivated
template<class SieveTraits, bool MT>
template<class LatP, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatP>)>
inline auto BlockOrthogonalSimHash<SieveTraits,MT>::transform_and_bitapprox_2nd_layer(LatP const &point)
    -> SimHashes<SieveTraits,MT>
{
  using std::abs;
  using std::max;
  unsigned int const dim = static_cast<unsigned int>(point.get_dim());

  using ET    = Get_AbsoluteCooType<LatP>;
  using Block = std::vector<ET>;  // will be of length dim.
  Block copy_of_point(dim);
  for (unsigned int i = 0; i < dim; ++i)
  {
    copy_of_point[i] = point.get_absolute_coo(i);
  }
  std::vector<Block> blocks(number_of_orthogonal_blocks, std::move(copy_of_point));

  for (uint_fast16_t i = 0; i < number_of_orthogonal_blocks; ++i)  // for each block:
  {
    for (uint_fast8_t j = 0; j < num_of_transforms; ++j)  // repeat num_of_transforms many times:
    {
      pmatrices[i][j].apply(blocks[i]);
      dmatrices[i][j].apply(blocks[i]);
      blocks[i] = fast_partial_walsh_hadamard(blocks[i], fast_walsh_hadamard_len);
    }
  }
  // put together the blocks into an array of bitsets.
  SimHashes<SieveTraits,MT> ret;
  for (unsigned int n = 0; n < sim_hash_num ; ++n)
  {
    // maximal entry of the outbut block.
    // Note: The walsh-hadamard transform used above is not normalized, so the elements are larger
    // than the original vector. That's why we do not compare to norm2() currently.
    // TODO : Reconsider the point where we change the sign. Currently, it's at
    // half the maximal entry. This may not be ideal.
    /*
    ET maxentry = 0;
    for (unsigned int m = 0; m < sim_hash_len; ++m)
    {
      unsigned int const flat_bit_count = n * sim_hash_len + m;
      maxentry = max(maxentry, static_cast<ET>(abs(blocks[flat_bit_count / dim][flat_bit_count % dim])));
    }
    */
    ET const norm2 = point.get_norm2();

    for (unsigned int m = 0; m < sim_hash_len; ++m)
    {
      // index of the bit currently considered
      // if we use only one level of indexing.
      unsigned int const flat_bit_count = n * sim_hash_len + m;
      ret[n][m] = ( blocks[flat_bit_count / dim][flat_bit_count % dim]
                  * blocks[flat_bit_count / dim][flat_bit_count % dim] * dim * 6 > norm2);
    }
//    std::cout << ret[0].count() << std::endl << std::flush;
  }
  return ret;
}
#endif // 0

/*
template<class SieveTraits, bool MT>
template<class LatP, TEMPL_RESTRICT_IMPL2(IsALatticePoint<LatP>)>
inline auto BlockOrthogonalSimHash<SieveTraits,MT>::transform_and_bitapprox_simple(LatP const &point)
    -> std::array< std::bitset<sim_hash_len>,num_of_levels >
{
  using ET    = Get_AbsoluteCooType<LatP>;
  ET res;

  std::array< std::bitset<sim_hash_len>, num_of_levels > ret;
  for (unsigned int i = 0; i < num_of_levels; ++i)
  {
    for (unsigned int j = 0; j < sim_hash_len; ++j)
    {
      for (unsigned int k=0; k < SimHash::num_of_coos; ++k)
      {
        uint_fast16_t coo = BlockOrthogonalSimHash<SieveTraits,MT>::rmatrices[i].get_value(i, k);
        res+=point[coo];
      }
      ret[i][j] = res>0;

    }

  }
  return ret;

}
 */


/** OLD UNUSED CODE
*/

/*
// helpers for the function below, do not use.
namespace GaussSieve{ namespace SimHash { namespace Helpers{
template<int SizeOfBitSet> struct MakeBitApprox_Helper
{
  static_assert(SizeOfBitSet >= 0, "Only for fixed-size bit-sets.");
  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  static inline std::bitset<SizeOfBitSet> compute_bitapproximation(LatP const &point)
  {
    auto const dim = point.get_dim();
    static_assert(dim == SizeOfBitSet, "size of bitset must equal fixed ambient dim of vector");
    std::bitset<SizeOfBitSet> ret;
    for (uint_fast16_t i = 0; i < dim; ++i)
    {
      ret[i] = (point.get_absolute_coo(i) >= 0) ? 1 : 0;
    }
    return ret;
  }
};

// variable-dim version:
template<> struct MakeBitApprox_Helper<-1>
{
  template<class LatP, TEMPL_RESTRICT_DECL2(IsALatticePoint<LatP>)>
  static inline boost::dynamic_bitset<> compute_bitapproximation(LatP const &point)
  {
    auto const dim = point.get_dim();
    boost::dynamic_bitset<> ret{ static_cast<size_t>(dim) };
//    ret.resize(dim);
    for (uint_fast16_t i = 0; i < dim; ++i)
    {
      ret[i] = (point.get_absolute_coo(i) >= 0) ? 1 : 0;
    }
    return ret;
  }
};
}}} // end namespaces GaussSieve::SimHash::Helpers
*/


/**
  compute_coordinate_wise_bitapproximation<length>(point)
  computes the coordinate wise 1-bit-approximation of the lattice point point.
  It returns a std::bitset<length> for lenght>= 0 or boost::dynamic_bitset for length == -1.
  Currently, the length parameter (if >=0) has to match the (fixed!) dimension of the lattice point.
  For length==1, the dynamic_bitset that is returned has length equal to the dimension of point.
*/

/*
namespace GaussSieve{ namespace SimHash{
template<int SizeOfBitSet, class LatP>
auto compute_coordinate_wise_bitapproximation(LatP const &point)
    -> decltype(  Helpers::MakeBitApprox_Helper<SizeOfBitSet>::compute_bitapproximation( std::declval<LatP>() )  )
{
  static_assert(IsALatticePoint<mystd::decay_t<LatP>>::value,"Not a lattice point.");
  return Helpers::MakeBitApprox_Helper<SizeOfBitSet>::compute_bitapproximation(point);
}
}}// end namespace GaussSieve::SimHash
*/


/*
namespace GaussSieve{ namespace SimHash{

// assume sim_hash2_len is a power-of-two.

// deprecated. Use fast_partial_walsh_hadamard instead.
template<class T>
[[deprecated]] inline std::vector<T> fast_walsh_hadamard_ext(std::vector<T> input, unsigned int len)
{
  std::vector<T> output (input.size());
  std::vector<T> tmp = input;
  uint_fast16_t i, j, s;

  for (i = (len >> 1); i > 0; i >>= 1)
  {
    for (j = 0; j < len; j++)
    {
      s = (j/i)%2;
      output[j]=input[(s?-i:0)+j]+(s?-1:1)*input[(s?0:i)+j];
    }
    //tmp = inp; inp = out; out = tmp;
    std::swap(input, output);
  }

  //lower matrix with 1's on the main diag
  for (i = len; i < output.size(); ++i)
  {
    output[i] = tmp[i];
  }

  return output;
}

}}  // end namespace GaussSieve::SimHash
*/
