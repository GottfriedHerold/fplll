//
//  FilteredPoint.h
//  
//
//  Created by Elena on 06/03/17.
//
//

#ifndef _FilteredPoint_h
#define _FilteredPoint_h

#include "sieve_common.h"
#include "LatticePoint2.h"

using namespace LatticeApproximations;

template <class ET> class FilteredPoint;

//template <class ET, bool insideMTList=false, int n_fixed=-1>
template <class ET>
class FilteredPoint
{
    public:
    
    FilteredPoint()=default;
    FilteredPoint(const FilteredPoint &Point) = default; // : NumVect<ET>::data(Point.data), norm2(Point.norm2) {}
    FilteredPoint(FilteredPoint &&Point) = default ;
    FilteredPoint(ApproxLatticePoint<ET> x, ET sc)
    {
        this->point = x;
        this->sc_prod = sc;
    }
    
    // if sc is int_32
    FilteredPoint(ApproxLatticePoint<ET> x, int32_t sc)
    {
        this->point = x;
        this->sc_prod = sc;
    }
    
    
    //FilteredPoint(ApproxLatticePoint x, ApproxLatticePoint p)
    
    
    FilteredPoint& operator=(FilteredPoint const &that) =default;
    FilteredPoint& operator=(FilteredPoint && that) =default;
    
    
    ~FilteredPoint() {}
    
    inline ApproxType * getApproxVector() const {return this->point.get_approx();}
    inline ET get_sc_prod() const {return sc_prod;}

    
private:
    //members
    ApproxLatticePoint<ET> point;
    ET sc_prod;
    
    
    
};

#endif
