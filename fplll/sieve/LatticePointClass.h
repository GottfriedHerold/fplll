#ifndef LATTICE_POINT_CLASS_H
#define LATTICE_POINT_CLASS_H

#include <type_traits>
#include "Utility.h"
#include <vector>

template<class ET,int nfixed>
class MyLatticePoint{
    
    using LatticePointType = true_type;
    using AuxDataType = Dimension<nfixed>;
    using ScalarProductReturnType = ET;
    
public:
    MyLatticePoint()=delete;
    MyLatticePoint(const MyLatticePoint &Point) = delete;
    MyLatticePoint(MyLatticePoint &&Point) = default ;
    
    MyLatticePoint& operator=(MyLatticePoint const &that) =default;
    MyLatticePoint& operator=(MyLatticePoint && that) =default;
    ~ExactLatticePoint() {}
    
    explicit MyLatticePoint(){
        
    };
    
    
    explicit MyLatticePoint(MatrixRow<ET> const & row) {
        data = (row.get_underlying_row()).get();
        update_norm2(data);
    };

    
    ET get_norm2() {return this->norm2;}
    
    
    
public:
    
    std::vector<ET> data;
    ET norm2;
};


template <class ET,int nfixed> MyLatticePoint<ET> add (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> MyLatticePoint<ET,nfixed> minus (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET, nfixed> const &B, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> MyLatticePoint<ET,nfixed> negate (MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> void scalar_mult (MyLatticePoint <ET,nfixed> &A, ET const & multiple, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> bool comapre_sc_product (MyLatticePoint<ET, fixed> const &A, MyLatticePoint<ET,nfixed> const &B,  ET target);
template <class ET,int nfixed> bool comapre_abs_sc_product (MyLatticePoint<ET, fixed> const &A, MyLatticePoint<ET,nfixed> const &B,  ET target);
template <class ET,int nfixed> ET compute_sc_product (MyLatticePoint<ET, fixed> const &A, MyLatticePoint<ET,nfixed> const &B);
template <class ET,int nfixed> MyLatticePoint<ET> make_copy (MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> auxdata);
#endif
