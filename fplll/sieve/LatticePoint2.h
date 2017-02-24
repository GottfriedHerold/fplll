//
//  LatticePoint.h -- new version
//
#ifndef LATTICE_VECTOR_CLASS2_H
#define LATTICE_VECTOR_CLASS2_H

#include "sieve_common.h" //needed (at least for convertions from MatrixRow (the header has to be revised);
#include "LatticePoint.h"
#include <type_traits>

//forward declarations:

template <class ET,bool insideMTList=false, int n_fixed=-1> //insideMTList: Is this point stored inside our Multithreaded list class? -> changes some data types, memory layout and ownership semantics.
                                                            //n : compile-time fixed dimension with n=-1 meaning dynamic.
class ApproxLatticePoint;

template<class ET, bool insideMTList, int n_fixed>
ostream & operator<< (ostream &os, ApproxLatticePoint<ET,insideMTList,n_fixed> const & appLP); //output.


namespace LatticeApproximations //helper types etc. enclosed in namespace
{

inline signed int get_exponent (Z_NR<mpz_t> const & val)     {return val!=0 ?val.exponent() : std::numeric_limits<signed int>::min();}     //returns smallest t s.t. abs(val) < 2^t (i.e. the bit-length for integral types)
inline signed int get_exponent (Z_NR<long> const & val)      {return val!=0 ?val.exponent() : std::numeric_limits<signed int>::min();}     //for val==0, we return the most negative value that fits into an int.
inline signed int get_exponent (Z_NR<double> const & val)    {return val!=0 ?val.exponent() : std::numeric_limits<signed int>::min();}     //Note that the function is already implemented in Z_NR< > - types.
template<class ET> [[ deprecated ("Using badly supported type") ]] signed int get_exponent (ET const & val);                               //Non-Z_NR< > - cases. Note that non-templates take precedence.

template<class ApproxType, class ET> ApproxType do_approximate( typename enable_if< is_same<ET, Z_NR<double> >::value, ET>::type const & val, signed int const delta);
template<class ApproxType, class ET> ApproxType do_approximate( typename enable_if< is_same<ET, Z_NR<mpz_t> >::value , ET>::type const & val, signed int const delta);
template<class ApproxType, class ET> ApproxType do_approximate( typename enable_if< is_same<ET, Z_NR<long> >::value , ET>::type const & val, signed int const delta);
// "if constexpr" (from C++17), you are badly needed...

using ApproxTypeNorm2 = int32_t; //determines bit-length of approximation, by making sure that n* (MAXAPPROX^2) must fit into it.
using ApproxType = int32_t; //at least half the size of the above. Same size might work better for vectorization.


//The above should be read as specialisations of
//template<class ET,class ApproxType>
//ApproxType do_approximate(ET const & val, signed int const delta );
//for E = Z_NR<...>

template<class ET>
class MaybeRational;

template<class ET> class MaybeRational {public: static bool constexpr val=true;};       //helper template that selects whether the underlying type *might* be able to represent values strictly between 0 and 1.
template<> class MaybeRational<Z_NR<long > >{public: static bool constexpr val=false;}; //this distinction just to serves to avoid correct, but needless approximations if the values are small enough to not need approximations in first place.
template<> class MaybeRational<Z_NR<mpz_t> >{public: static bool constexpr val=false;}; //(otherwise, we would pad with zeros from the right(!), which is correct (but useless and hurts readability for debug outputs).
template<> class MaybeRational<Z_NR<double> >{public: static bool constexpr val=true;};
//template<class ET>

inline ApproxTypeNorm2 compute_sc_prod(ApproxType const * const arg1, ApproxType const * const arg2, unsigned int len);
template<class ET>
inline bool Compare_Sc_Prod(ApproxLatticePoint<ET,false,-1> const & arg1, ApproxLatticePoint<ET,false,-1> const & arg2, ApproxTypeNorm2 abslimit, int limit_exp, int dim);
}




//end of forward declarations.

//TODO: Consider Renaming AppproxLatticePoint -> LatticePoint and
//                        LatticePoint        -> LatticePointDetails
//to better match the semantics.

template <class ET>
class ApproxLatticePoint<ET, false, -1>     //Approx. Lattice Point. Stores approximations to a lattice point LP in the form LP = 2^length_exponent * (*approx).
                                            //i.e. essentially like a floating-point type, but with an exponent that is shared between all coordinates.
                                            //The underlying exact value can be obtained by get_details.

                            //Note that the pointer/handle to the underlying LatticePoint is owning.
                            //The reasong for this strange design choice (rather doing it the other way round, i.e. attaching the approx. to the real thing)
                            //is that the memory layout of the approx. data is completly under our control, does not depend on <ET>, does not grow and
                            //every potential sequence of bytes represents at least some valid object.
                            //In a multi-threaded environment, the latter allows us to forego some locking / relax memory order requirements at the expense of having approximations that actually are wrong.
                            //Dealing with this may be cheaper than having memory fences at every read. Of course, reading from the "true" values requires some more care in that case.

{
    public:
    using ApproxType        = LatticeApproximations::ApproxType;
    using ApproxTypeNorm2   = LatticeApproximations::ApproxTypeNorm2;
    using DetailType        = LatticePoint<ET>;

    public: //consider making some constructors private and befriend the list class(es).
    ApproxLatticePoint() : length_exponent(0),approx(nullptr), approxn2(0), details(nullptr) {}; //should only ever be used in move semantics
    //explicit ApproxLatticePoint(int n) : //n=ambient dimension = length of vector
    //    length_exponent(0), approx(new ApproxType[n]),approxn2(0),details(new DetailType (n) ) {update_approximation();} ;
    ApproxLatticePoint(DetailType const & LP): //creates an approx. LP structure holding a copy of LP.
        length_exponent(0),approx(new ApproxType[LP.size()] ),approxn2(0),details(new DetailType (LP)) {update_approximation();};
    //explicit ApproxLatticePoint(DetailType const * const LPp): //creates an approx LP structure that gains ownership of *LPp. Use only on dynamically allocated LPs.
    //    length_exponent(0),approx(new ApproxType[LPp->size()]),approxn2(0),details(LPp) {update_approximation();};

    ApproxLatticePoint(ApproxLatticePoint const &other):
        length_exponent(other.length_exponent),approx(nullptr),approxn2(other.approxn2),details(nullptr)
    {
        if (other.details!=nullptr)
        {
            int n = (other.details)->size();
            approx= new ApproxType[n];
            for(int i=0;i<n;++i) approx[i]=other.approx[i];
            details = new DetailType (*other.details);
        }
    }
    ApproxLatticePoint(ApproxLatticePoint && other) :
        length_exponent(other.length_exponent), approx(other.approx), approxn2(other.approxn2),details(other.details)  {other.invalidate();}; //invalidation should not be neccessary.
    //ApproxLatticePoint & operator= (ApproxLatticePoint const & other) = delete;
    ApproxLatticePoint & operator= (ApproxLatticePoint other)
    {
        length_exponent = other.length_exponent;
        swap(approx,other.approx);
        approxn2=other.approxn2;
        swap(details,other.details);
        return *this;
    }

    ~ApproxLatticePoint(){delete details; delete approx;};
    ApproxType * get_approx()           const               {return approx;};
    ApproxType   get_approx_norm2()     const               {return approxn2;};
    signed int get_length_exponent()  const                 {return length_exponent;};
    DetailType get_details()    const                       {return *details;};
    //LatticePoint<ET> get_details()    const                 {return *details;};
    DetailType & get_details_ref()                          {return *details;}; //technically, this is const, as it does not change the pointer.
    DetailType const * get_details_ptr()  const             {return details;};
    DetailType const * get_details_ptr_r()const             {return details;};
    DetailType * get_details_ptr_rw() const                 {return details;};
    unsigned int get_dim() const                            {return details->size();};
    bool operator< (ApproxLatticePoint const &other ) const {return (this->details->norm2 < other.get_details_ptr()->norm2);};
    private: //internal member functions
    void update_approximation(); //computes approximation from *approx. Assumes memory is already allocated.
    void invalidate(){approx= nullptr;details=nullptr;}
    //friend
    private: //internal data
    signed int length_exponent; //note : May be negative
    ApproxType *approx;
    ApproxTypeNorm2 approxn2;
    //The approximation itself is given by 2^length_exponent * approx
    //The approximation to the norm^2 is given by 2^2length_exponent * approxn2
    //Note that we need to care about overflows here by truncating accordingly.
    DetailType *details; //actual lattice point structure.
};

//TODO: Implementation will be changed
template<class ET, bool insideMTList, int n_fixed>
ostream & operator<< (ostream &os, ApproxLatticePoint<ET,insideMTList,n_fixed> const & appLP)
{
    os << appLP.get_details() << endl;
    os << "Approximation:[";
    for (unsigned int i=0;i<appLP.get_dim(); ++i)
    {
        os<<(appLP.get_approx())[i]<<","; //gives trailing comma after last arg. Too lazy to fix.
    }
    os<<"]";
    os<<" with multiplier 2^" << appLP.get_length_exponent() << endl;
    os<<"norm2-Approx = " << appLP.get_approx_norm2() << endl;
    return os;
}


template <class ET>
void ApproxLatticePoint<ET,false, -1>::update_approximation()
{
    assert(details!=nullptr);
    int const n = details->size(); //should be unsigned
    assert(n!=0);
    static_assert (std::numeric_limits<ApproxType>::is_specialized, "bad ApproxType");
    static_assert (std::numeric_limits<ApproxTypeNorm2>::is_specialized,"bad ApproxTypeNorm2");
    static_assert (std::numeric_limits<ApproxTypeNorm2>::digits /2 <= std::numeric_limits<ApproxType>::digits, "bad Types"); //ensures that we have enough bits in ApproxType. Note that /2 rounds down, which is correct.
    static_assert ( (! is_same<ET,Z_NR<mpz_t>>::value) || ( ( std::numeric_limits<ApproxTypeNorm2>::digits / 2 ) <= std::numeric_limits<long>::digits )," " ); //otherwise, mpz_get_si in the GMP class may not work correctly.
    unsigned int const max_bits = ( std::numeric_limits<ApproxTypeNorm2>::digits - floor(log2(n)) )/2;
    //may use __builtin_clz on GCC for improved efficiency, floor(log2(n)) is ridiculous...
    signed int number_length = std::numeric_limits<signed int>::min();
    for(int i=0; i<n;++i)
    {
        number_length = std::max(number_length, LatticeApproximations::get_exponent ( (*details)[i] ) );
    }
    if(number_length == std::numeric_limits<signed int>::min()) // *details is all-zero vector. should never happen in the algorithm.
    {
        cerr << "Warning: approximating all-zero vector";
        length_exponent=0;
        for(int i=0;i<n;++i){approx[i]=0;}
        approxn2=0;
        return;
    }
    else      // length_exponent is minimal, s.t. abs ( (*details)[i] / 2^length_exponent ) < 1 for all i
    {
        for(int i=0;i<n;++i)
        {
            length_exponent = number_length - max_bits;
            if(LatticeApproximations::MaybeRational<ET>::val == false) {length_exponent = max(0,length_exponent);} //constexpr if, actually...
            approx[i] = LatticeApproximations::do_approximate<ApproxType,ET> ( (*details)[i], length_exponent );
        }
        approxn2 = LatticeApproximations::do_approximate<ApproxTypeNorm2,ET> ( (*details).get_norm2(),2*length_exponent   );
    }
}

template<class ApproxType, class ET>
ApproxType LatticeApproximations::do_approximate( typename enable_if< is_same<ET, Z_NR<mpz_t> >::value , ET>::type const & val, signed int const delta)
//ApproxType do_approximate(Z_NR<mpz_t> const & val, signed int const delta )
{
    assert(delta>=0);
    Z_NR<mpz_t> temp;
    temp.div_2si (val,delta);
    return temp.get_si();
}

//Note that nr_Z_l 's div_2si uses >> and << operators on signed values.
//This would lead to undefined behaviour (as explicitly stated in nr_Z) in our case as of C++17 (although it would work with most compilers).
//(This is the single reason we don't default to Z_NR's capabilities here.


template<class ApproxType, class ET>
ApproxType LatticeApproximations::do_approximate( typename enable_if< is_same<ET, Z_NR<long> >::value , ET>::type const & val, signed int const delta)
//ApproxType do_approximate(Z_NR<long> const & val, signed int const delta)
{
    assert(delta>=0);
    return ( val.get_data() / (2L << delta ) );
}

template<class ApproxType, class ET>
ApproxType LatticeApproximations::do_approximate( typename enable_if< is_same<ET, Z_NR<double> >::value , ET>::type const & val, signed int const delta)
//Read : ApproxType do_approximate(Z_NR<double> const & val, signed int const delta)
{
    return (std::ldexp(val.get_data(),-delta) ); //Note : Conversion from double to integral type rounds towards zero. This is needed to prevent overflows.
}


/*get_exponent :
returns minimal r, s.t. |val| < 2^r (or INTMIN, if val = 0)
*/

template <class ET> //fallback version, should never be called anyway.
[[ deprecated ("Using badly supported type") ]] signed int LatticeApproximations::get_exponent(ET const & val)
{

    if(val==0) return std::numeric_limits<signed int>::min();
    ET absval = abs(val);
    if(absval >= 1)
    {
        signed int res=0;
        while(absval>=1)
        {
            ++res;
            absval = absval /2;
        }
        return res;
    }
    else
    {
        signed int res=1;
        while(absval<1)
        {
            --res;
            absval = absval*2;
        }
        return res;
    }
}

inline LatticeApproximations::ApproxTypeNorm2 LatticeApproximations::compute_sc_prod(ApproxType const * const arg1, ApproxType const * const arg2, unsigned int len)
{
    ApproxTypeNorm2 res=0;
    for(unsigned int i=0;i<len;++i)
    {
        res+=arg1[i]*arg2[i];
    }
    return res;
}

template<class ET>
inline bool LatticeApproximations::Compare_Sc_Prod(ApproxLatticePoint<ET,false,-1> const & arg1, ApproxLatticePoint<ET,false,-1> const & arg2, ApproxTypeNorm2 abslimit, int limit_exp, int dim)
{
    int rel = arg1.get_length_exponent() + arg2.get_length_exponent() - limit_exp;
    ApproxTypeNorm2 sc = abs(compute_sc_prod(arg1.get_approx(),arg2.get_approx(),dim));
    if(rel > 0)
    {
        return sc > (abslimit >> rel);
    }
    else
    {
        return (sc >> -rel) > abslimit;
    }
}


#endif // LATTICE_VECTOR_CLASS_2H
