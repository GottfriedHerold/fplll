#ifndef COMPRESSED_LATTICE_POINT_H
#define COMPRESSED_LATTICE_POINT_H

template<class ET, bool MT, int nfixed>
class CompressedPoint
{
    public:
    using DetailType = ExactLatticePoint<ET,nfixed>;
    CompressedPoint() : approx_data(), details(nullptr) {};                                                                 //creates a dummy, invalid point
    CompressedPoint(CompressedPoint && old) : approx_data(std::move(old.approx_data)), details(old.details) {old.details = nullptr;}; //move semantics
    CompressedPoint(CompressedPoint const & old) = delete; //no copying
    explicit CompressedPoint(ExactLatticePoint<ET,nfixed> & exact_point);     //creates a CP from an exact point
    //explicit CompressedPoint(ExactLatticePoint && exact_point) = delete;        //creates a CP from an exact point, move semantics
    ~CompressedPoint();

    ApproximateLatticePoint<ET,nfixed>  const & access_approximation_r() const {return approx_data;};
    typename ApproximateLatticePoint<ET,nfixed>::ApproxTypeNorm2 get_approx_norm2_mantissa() const {return approx_data.get_norm2_mantissa();};
    signed int get_length_exponent() const {return approx_data.get_length_exponent();};
    DetailType get_exact_point() const;                 //returns a copy of the underlying exact point

    protected:
    ApproximateLatticePoint<ET,nfixed> approx_data;
    DetailType * details;
};


template<class ET, bool MT, int nfixed> CompressedPoint<ET,MT,nfixed>::CompressedPoint(ExactLatticePoint<ET,nfixed> & exact_point)
:approx_data(exact_point), details( &exact_point)
{
}


template<class ET, bool MT, int nfixed> typename CompressedPoint<ET,MT,nfixed>::DetailType CompressedPoint<ET,MT,nfixed>::get_exact_point() const
{
    return * details;
}

template<class ET, bool MT, int nfixed> CompressedPoint<ET,MT,nfixed>::~CompressedPoint()
{
    delete details;
}
#include "CompressedPoint.cpp"
#endif
