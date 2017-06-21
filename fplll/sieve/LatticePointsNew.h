/* New header file for classes storing lattice points. */

#ifndef LATTICE_VECTOR_CLASS_NEW_H
#define LATTICE_VECTOR_CLASS_NEW_H

#include "sieve_common.h" //needed (at least for convertions from MatrixRow (the header has to be revised);
#include <string>
#include <numeric>

//forward declarations:


template <class ET,int nfixed=-1> class ExactLatticePoint;      //stores an (exact) lattice point, together with its (squared) length.
                                                                //nfixed denotes a fixed dimension (at compile-time) of vectors, with -1 meaning runtime-determined.
                                                                //not currently implemented.

//TODO: Move to User.
//template <class ET,int nfixed=-1> class IsLongerExactVector_class;  //class wrapper to compare vectors by length.
                                                                    //Needs to be wrapped in a class to work seamlessly with some STL containers.
                                                                    //Alternatively, we could use lambdas...
template <class ET,int nfixed=-1> void compute_exact_sc_product (ET &result, const ExactLatticePoint<ET,nfixed> &p1, const ExactLatticePoint<ET,nfixed> &p2);   //computes the scalar product of two exact lattice points.
template <class ET,int nfixed=-1> ET exact_scalar_product(ExactLatticePoint<ET,nfixed> const &p1, ExactLatticePoint<ET,nfixed> const & p2);                     //same, but returns result instead of storing in first arg.
template <class ET,int nfixed=-1> void scalar_mult (ExactLatticePoint<ET,nfixed> &A, ET const & multiple); //A = A*multiple
template <class ET,int nfixed=-1> std::ostream & operator<<(std::ostream &os, ExactLatticePoint<ET,nfixed> const & exact_point); //printing
//template<class ET> istream & operator>>(istream &is, ExactLatticePoint & exact_point); //reading (may also be used by constructor from istream)
//FIXME: operator >> is NOT IMPLEMENTED YET

//FIXME: Enable or use a constructor
//template<class ET,int nfixed=-1> ExactLatticePoint<ET,nfixed> conv_matrixrow_to_lattice_point (MatrixRow<ET> const &row);

//-----------------------------------

/*
    ApproximateLatticePoint stores an approximation to a lattice point in the form
    2^exponent * (mantissa), where exponent is shared among coordinates and
    (mantissa) is a vector of type LatticeApproximations::ApproxType (should be approx 16-32 bits)
    We also store a copy of the norm^2 of the mantissa (so the real norm is 2^2exponent * norm_of_mantissa
    Note that is does not really depend much on the underlying type ET, which is only used to select the correct conversion routines and to know whether exponent may be negative.
*/

template<class ET, int nfixed=-1> class ApproximateLatticePoint;      //This class represents an approximation to a lattice point.
                                                                        //Note that it is detached from whatever it approximates.
template<class ET, int nfixed=-1> std::ostream& operator<< (std::ostream &os, ApproximateLatticePoint<ET,nfixed> const & approx_point); //output.
//TODO: Input

/* Creates a stand-alone Approximation from an exact point */
template<class ET, int nfixed> ApproximateLatticePoint<ET,nfixed> create_detached_approximation(ExactLatticePoint<ET,nfixed> const & exact_point);

//----------------------------------------//

//This class represents an approximate point together with some extra data (called Details) that allows to recover the exact point again.
//Note that access to the exact point may be slow (in multi-threading, it might involve locks or fences), so working with the approximation is prefered.
//IMPORTANT:


template<class ET, bool MT, int nfixed=-1> class CompressedPoint;




namespace LatticeApproximationsNew //internal helper types for approximations etc. enclosed in namespace
{
inline signed int get_exponent (Z_NR<mpz_t> const & val)     {return val!=0 ?val.exponent() : std::numeric_limits<signed int>::min();}     //returns smallest t s.t. abs(val) < 2^t (i.e. the bit-length for integral types)
inline signed int get_exponent (Z_NR<long> const & val)      {return val!=0 ?val.exponent() : std::numeric_limits<signed int>::min();}     //for val==0, we return the most negative value that fits into an int.
inline signed int get_exponent (Z_NR<double> const & val)    {return val!=0 ?val.exponent() : std::numeric_limits<signed int>::min();}     //Note that the function is already implemented in Z_NR< > - types.
template<class ET> [[ deprecated ("Using badly supported type") ]] signed int get_exponent (ET const & val);                               //Non-Z_NR< > - cases. Note that non-templates take precedence.


//ApproxType do_approximate(X val, int delta) returns 2^(-delta)*val, converted to appropriate ApproxType.
//Implementation depends on X

template<class ApproxType, class ET> ApproxType do_approximate( typename enable_if< is_same<ET, Z_NR<double> >::value, ET>::type const & val, signed int const delta);
template<class ApproxType, class ET> ApproxType do_approximate( typename enable_if< is_same<ET, Z_NR<mpz_t > >::value, ET>::type const & val, signed int const delta);
template<class ApproxType, class ET> ApproxType do_approximate( typename enable_if< is_same<ET, Z_NR<long  > >::value, ET>::type const & val, signed int const delta);

using ApproximationNorm2Type = int32_t; //determines bit-length of approximation, by making sure that n* (MAXAPPROX^2) must fit into it.
using ApproximationEntriesType = int32_t; //at least half the size of the above. Same size might work better for vectorization.

template<class ET> class MaybeRational; //Compile-time function of ET that tells whether ET might represent values between 0 and 1
template<class ET> class MaybeRational       {public: static bool constexpr val=true;};  //helper template that selects whether the underlying type *might* be able to represent values strictly between 0 and 1.
template<> class MaybeRational<Z_NR<long > > {public: static bool constexpr val=false;}; //this distinction just to serves to avoid correct, but needless approximations if the values are small enough to not need approximations in the first place.
template<> class MaybeRational<Z_NR<mpz_t> > {public: static bool constexpr val=false;}; //(otherwise, we would pad with zeros from the right(!), which is correct (but useless and hurts readability for debug outputs).
template<> class MaybeRational<Z_NR<double> >{public: static bool constexpr val=true;};
//template<class ET>

//inline ApproxTypeNorm2 compute_sc_prod(ApproxType const * const arg1, ApproxType const * const arg2, unsigned int len);

//template<class ET>
//inline bool Compare_Sc_Prod(ApproxLatticePoint<ET,false,-1> const & arg1, ApproxLatticePoint<ET,false,-1> const & arg2, ApproxTypeNorm2 abslimit, int limit_exp, int dim);
};

namespace GaussSieve
{
//computes scalar products: mantissa variant only computes the scalar product of the mantissas (i.e. it lacks an 2^{exponent1+exponent2} - factor)
//                          vec version computes scalar products between two arrays.
template<class ET,int nfixed>   LatticeApproximationsNew::ApproximationNorm2Type compute_mantissa_sc_product(ApproximateLatticePoint<ET,nfixed> const & arg1, ApproximateLatticePoint<ET,nfixed> const & arg2, int const dim);
                                LatticeApproximationsNew::ApproximationNorm2Type compute_vec_sc_products(LatticeApproximationsNew::ApproximationEntriesType const * const arg1, LatticeApproximationsNew::ApproximationEntriesType const * const arg2, int const dim);
//FIXME: templating by nfixed does not work here. Correct way is to make ApproxEntryType* a template.
//The error message is weird and makes me suspect compiler bugs (the issue is with template function signature not depending on template param, but it give namespace errors (GCC)
};



/*actual class declarations in separate files*/

#include "ExactLatticePoint.h"
#include "ApproximateLatticePoint.h"
#include "CompressedPoint.h"

LatticeApproximationsNew::ApproximationNorm2Type GaussSieve::compute_vec_sc_products(LatticeApproximationsNew::ApproximationEntriesType const * const arg1, LatticeApproximationsNew::ApproximationEntriesType const * const arg2, int const dim)
{
    return std::inner_product(arg1, arg1+dim, arg2,0);
}

//template<int nfixed> LatticeApproximationsNew::ApproximationNorm2Type GaussSieve::compute_vec_sc_products<-1>(LatticeApproximationsNew::ApproximationEntriesType const * const arg1, LatticeApproximationsNew::ApproximationEntriesType const * const arg2, int const dim)
//{
//    return std::inner_product(arg1, arg1+nfixed,arg2,0);
//}

template<class ET,int nfixed>  LatticeApproximationsNew::ApproximationNorm2Type GaussSieve::compute_mantissa_sc_product(ApproximateLatticePoint<ET,nfixed> const & arg1, ApproximateLatticePoint<ET,nfixed> const & arg2, int const dim)
{
    return GaussSieve::compute_vec_sc_products(arg1.access_vectors_mantissa(),arg2.access_vectors_mantissa(),dim);
}

#include "LatticePointsNew.cpp"

//instantiate for appropriate types
template class ExactLatticePoint<Z_NR<long>>;
template class ExactLatticePoint<Z_NR<double>>;
template class ExactLatticePoint<Z_NR<mpz_t>>;
template class ApproximateLatticePoint<Z_NR<long>>;
template class ApproximateLatticePoint<Z_NR<double>>;
template class ApproximateLatticePoint<Z_NR<mpz_t>>;
template class CompressedPoint<Z_NR<long>,false>;
template class CompressedPoint<Z_NR<mpz_t>,false>;
template class CompressedPoint<Z_NR<double>,false>;

#endif
