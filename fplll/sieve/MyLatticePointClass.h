#ifndef MY_LATTICE_POINT_CLASS_H
#define MY_LATTICE_POINT_CLASS_H

#include <iostream>
#include <type_traits>
#include "Utility.h"
#include <vector>

template<class ET,int nfixed> class MyLatticePoint;
//template <class ET,int nfixed> ET compute_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata);




template<class ET,int nfixed>
class MyLatticePoint{
    
    using LatticePointType = true_type;
    using AuxDataType = Dimension<nfixed>;
    using ScalarProductReturnType = ET;
    
public:
    MyLatticePoint()=delete;
    MyLatticePoint(MyLatticePoint const &Point) = delete;
    MyLatticePoint(MyLatticePoint &&Point) = default ;
    
    MyLatticePoint& operator=(MyLatticePoint const &that) =delete;
    MyLatticePoint& operator=(MyLatticePoint && that) =default;
    ~MyLatticePoint() {}
    
    
    //explicit MyLatticePoint(Dimension<nfixed> dim={}, IgnoreArg<AuxDataType const &> auxdata = {}){ 'auxdata = {}' <- ERRS!
    explicit MyLatticePoint(Dimension<nfixed> dim, IgnoreArg<AuxDataType const &> auxdata ){
        data = std::vector<ET>(dim.dim);
        ET norm;
        norm = 0;
        norm2 = norm;
    };
    
    
    explicit MyLatticePoint(MatrixRow<ET> const & row, Dimension<nfixed> const & dim) {
        data = (row.get_underlying_row()).get();
        update_norm2( dim);
    };
    
    
    //In order to be able to make a copy
    MyLatticePoint(MyLatticePoint<ET, nfixed> const & point_, Dimension<nfixed> dim)
    {
        int n = dim.dim;
        data =  std::vector<ET>(n);
        //std::memcpy(&data, &point_.data, n*sizeof(ET));
        norm2 = point_.norm2;
        data = point_.data;
        //norm2 =point_.norm2;
    }
    
    
    void update_norm2(Dimension<nfixed> const & dim)
    {
        this->norm2 = compute_sc_product(*this, *this, dim);
    };
    
    ET get_norm2() {return this->norm2;};
    
    //friend std::ostream & operator<< <ET, nfixed> (std::ostream &os, MyLatticePoint<ET,nfixed> const &A);
    
    void print (std::ostream &os, Dimension<nfixed> const & dim)
    {
        os << "[";
        int n = dim.dim;
        for (int i =0; i<n; ++i)
        {
            os << this->data[i] << " " ;
        }
        
        os <<"]" << endl;
    }
    
    
    
public:
    std::vector<ET> data;
    ET norm2;
};


template <class ET,int nfixed> MyLatticePoint<ET, nfixed> add (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> MyLatticePoint<ET,nfixed> sub (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET, nfixed> const &B, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> MyLatticePoint<ET,nfixed> negateP (MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata);


template <class ET,int nfixed> MyLatticePoint<ET,nfixed> scalar_mult (MyLatticePoint <ET,nfixed> &A, ET const & multiple, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> bool comapre_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B,  ET const & target);
template <class ET,int nfixed> bool comapre_abs_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B,  ET const & target);
template <class ET,int nfixed> ET compute_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> MyLatticePoint<ET, nfixed> make_copy (MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata);
//template <class ET,int nfixed> MyLatticePoint<ET, nfixed> void print (std::ostream &os = cout, MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata);


#endif
