#ifndef COMPRESSED_LATTICE_POINT_H
#define COMPRESSED_LATTICE_POINT_H

template<class ET, bool MT, int nfixed>
class CompressedPoint
{
    public:
    using DetailType = ExactLatticePoint<ET,nfixed> *;
    CompressedPoint() : approx_data(), details(nullptr) {};
    CompressedPoint(CompressedPoint && old) : approx_data(old.approx_data), detail(old.details) {old.details = nullptr;}; //move semantics
    CompressedPoint(CompressedPoint const & old) = delete; //no copying
    explicit CompressedPoint(ExactLatticePoint const & exact_point);    //creates a CP from an exact point
    explicit CompressedPoint(ExactLatticePoint && exact_point);         //creates a CP from an exact point, move semantics



    ApproximateLatticePoint<ET,nfixed>  get_approximation() const {return approx_data;};
    typename ApproximateLatticePoint<ET,nfixed>::ApproxTypeNorm2 get_approx_norm2_mantissa() const {return approx_data.get_norm2_mantissa();};
    signed int get_length_exponent() const {return approx_data.get_length_exponent();};


    protected:
    ApproximateLatticePoint<ET,nfixed> approx_data;
    DetailType details;
}


template<class ET, bool MT, int nfixed> CompressedPoint<ET,MT,nfixed>::CompressedPoint(ExactLatticePoint const & exact_point)
template<class ET, bool MT, int nfixed> CompressedPoint<ET,MT,nfixed>::CompressedPoint(ExactLatticePoint && exact_point)
CompressedPoint(ExactLatticePoint const & exact_point);    //creates a CP from an exact point
CompressedPoint(ExactLatticePoint && exact_point);         //creates a CP from an exact point, move semantics
#include "CompressedPoint.cpp"
#endif


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
                            //Even in single-threaded environment, this potentially allows better custom allocation.

{

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
    bool operator<=(ApproxLatticePoint const &other ) const {return (this->details->norm2 <=other.get_details_ptr()->norm2);};
    bool operator> (ApproxLatticePoint const &other ) const {return (this->details->norm2 > other.get_details_ptr()->norm2);};
    bool operator>=(ApproxLatticePoint const &other ) const {return (this->details->norm2 >=other.get_details_ptr()->norm2);};
    //private: //internal member functions
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
