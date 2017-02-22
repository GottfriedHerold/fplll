//
//  LatticePoint.h -- new version
//
#ifndef LATTICE_VECTOR_CLASS2_H
#define LATTICE_VECTOR_CLASS2_H

#include "sieve_common.h" //needed (at least for convertions from MatrixRow (the header has to be revised);
#include "LatticePoint.h"

//forward declarations:

template <class ET,bool insideMTList, int n_fixed=-1> //insideMTList: Is this point stored inside our Multithreaded list class? -> changes some data types, memory layout and ownership semantics.
                                                      //n : compile-time fixed dimension with n=-1 meaning dynamic.
class ApproxLatticePoint;

template <class ET>
class ApproxLatticePoint<ET, false, -1>    //Approx. Lattice Point. Note that the pointer/handle to the underlying LatticePoint is owning.
                            //The reasong for this strange design choice (rather doing it the other way round, i.e. attaching the approx. to the real thing)
                            //is that the memory layout of the approx. data is completly under our control, does not depend on <ET>, does not grow and
                            //every potential sequence of bytes represents at least some valid object.
                            //In a multi-threaded environment, the latter allows us to forego some locking / relax memory order requirements at the expense of having approximations that actually are wrong.
                            //Dealing with this may be cheaper than having memory fences at every read. Of course, reading from the "true" values requires some more care in that case.

{
    public:
    using ApproxType = uint32_t;
    using ApproxTypeNorm2 = uint32_t;

    public: //consider making some constructors private and befriend the list class(es).
    ApproxLatticePoint() : length_exponent(0),approx(nullptr), approxn2(0), details(nullptr) {}; //should only ever be used in move semantics
    ApproxLatticePoint(int n) : length_exponent(0), approx(new ApproxType[n]),approxn2(0),details(new LatticePoint<ET> (n) ) {update_approximation();} ; //n=ambient dimension = length of vector
    ApproxLatticePoint(ApproxLatticePoint const &other)=delete;
    ApproxLatticePoint(ApproxLatticePoint && other)=default;
    ApproxLatticePoint & operator= (ApproxLatticePoint const & other) = delete;
    ApproxLatticePoint & operator= (ApproxLatticePoint &&other) = default;
    ~ApproxLatticePoint(){delete details; delete approx;};
    ApproxType * get_approx()           const               {return approx;};
    ApproxType   get_approx_norm2()     const               {return approxn2;};
    unsigned int get_length_exponent()  const               {return length_exponent;};
    LatticePoint<ET> get_details_r()    const               {return *details;};
    LatticePoint<ET> & get_details_rw()                     {return *details;}; //technically, this is const, as it does not change the pointer.
    LatticePoint<ET> const * get_details_ptr_r()const       {return details;};
    LatticePoint<ET> * get_details_ptr_rw() const           {return details;};

    private: //internal member functions
    void update_approximation(); //computes approximation from *approx. Assumes memory is already allocated.

    private: //internal data
    unsigned int length_exponent;
    ApproxType *approx;
    ApproxTypeNorm2 approxn2;
    //The approximation itself is given by 2^length_exponent * approx
    //The approximation to the norm^2 is given by 2^2length_exponent * approxn2
    //Note that we need to care about overflows here by truncating accordingly.
    LatticePoint<ET> *details; //actual lattice point structure.
};



#endif // LATTICE_VECTOR_CLASS_2H
