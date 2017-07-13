#ifndef MY_LATTICE_POINT_CLASS_H
#define MY_LATTICE_POINT_CLASS_H

#include <iostream>
#include <type_traits>
#include "Utility.h"
#include <vector>
#include "LatticePointConcept.h"


template<class ET,int nfixed> class MyLatticePoint;

template<class ET,int nfixed>
class ImplementationTraits<MyLatticePoint<ET,nfixed> > : public ImplementationTraitsBase
{
    public:
    using AuxDataType = Dimension<nfixed>;
    using ScalarProductReturnType = ET;
};

//template<class ET,int nfixed> class MyLatticePoint : public GeneralLatticePoint< MyLatticePoint<ET, nfixed>  >
//template <class ET,int nfixed> ET compute_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata);


template<class ET,int nfixed>
class MyLatticePoint : public GeneralLatticePoint< MyLatticePoint<ET, nfixed> >
{
    public:
    using LatticePointTag = true_type;
    using AuxDataType = typename ImplementationTraits<MyLatticePoint>::AuxDataType;
    using ScalarProductReturnType = ET;

public:
    MyLatticePoint()=delete;
    MyLatticePoint(MyLatticePoint const &point) = delete;
    MyLatticePoint(MyLatticePoint &&point) = default;
    MyLatticePoint& operator=(MyLatticePoint const & other) = delete;
    MyLatticePoint& operator=(MyLatticePoint && other) = default;

    ~MyLatticePoint() {}

    explicit MyLatticePoint(Dimension<nfixed> dim={}, IgnoreArg<AuxDataType const &> auxdata = {})
    {
    //{ 'auxdata = {}' <- ERRS!
    //explicit MyLatticePoint(Dimension<nfixed> dim, IgnoreArg<AuxDataType const &> auxdata ){
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
    explicit MyLatticePoint(MyLatticePoint<ET, nfixed> const & point_, Dimension<nfixed> dim)
    {
        data =  std::vector<ET>(dim);
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

    std::ostream & write_to_stream (std::ostream &os, Dimension<nfixed> const & dim)
    {
        os << "[";
        for (unsigned int i =0; i<dim; ++i)
        {
            os << this->data[i] << " " ;
        }

        os <<"]" << endl;
        return os;
    }

    //friend std::ostream & operator<< <ET, nfixed> (std::ostream &os, MyLatticePoint<ET,nfixed> const &A);

    MyLatticePoint make_copy (Dimension<nfixed> const & auxdata={}) const;
    
    ScalarProductReturnType get_norm2 (AuxDataType const & aux_data)
    {
        return norm2;
    }


public:
    std::vector<ET> data;
    ET norm2;
};




template <class ET,int nfixed> MyLatticePoint<ET, nfixed> add (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> MyLatticePoint<ET,nfixed> sub (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET, nfixed> const &B, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> MyLatticePoint<ET,nfixed> negateP (MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata);


template <class ET,int nfixed> MyLatticePoint<ET,nfixed> scalar_mult (MyLatticePoint <ET,nfixed> &A, ET const & multiple, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> bool compare_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, ET const & target,  Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> bool compare_abs_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, ET const & target,  Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> ET compute_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata);
//template <class ET,int nfixed> MyLatticePoint<ET, nfixed> void print (std::ostream &os = cout, MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata);


#endif
