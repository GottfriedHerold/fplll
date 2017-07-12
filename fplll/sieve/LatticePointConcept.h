#ifndef LATTICE_POINT_CONCEPT_H
#define LATTICE_POINT_CONCEPT_H

//TODO: Remove the concept idea below, we replace it by CRTP (i.e. inheritance without virtual functions, using templates)

#include "Utility.h"


/**

Lattice points represent points in the lattice.

We treat LatticePoints as a concept, in the sense that a lattice point is any class that implements the following interface:
(or a subset thereof)
(This is WIP, not sure what should be mandatory)
//Note that functions might also have a different interface with arguments that are convertible from AuxDataType


class LatticePoint{
public: //typedefs
using LatticePoint_Tag = true_type //required, because we static_assert that this is there or use it for SFINAE
using AuxDataType = ... //optional public typedef that specifies some class, defaults to IgnoreAnyArg
using ScalarProductReturnType //type returned by scalar products and norms

public: //member functions:
LatticePoint make_copy(AuxDataType& const aux_data); //actually makes a copy
LatticePoint(LatticePoint &&); //move constructor
LatticePoint(LatticePoint const &)=delete //copy constructor deleted. Always make copies explicitly
LatticePoint& operator=(LatticePoint const &old)=delete;
LatticePoint& operator=(LatticePoint && old);
ScalarProductReturnType get_norm2(AuxDataType const &aux_data);
void negate(AuxDataType); //negate function.
explicit LatticePoint(Dimension<nfixed> dim, AuxDataType const& aux_data); //creates an unitialized point of dimension dim.
void make_zero(AuxDataType);


}

//non-member functions:
LatticePoint add(LatticePoint const & A, LatticePoint const & B, A::AuxDataType);
LatticePoint subtract(LatticePoint const &A, LatticePoint const & B, A::AuxDataType);
A::ScalarProductReturnType compute_sc_product(LatticePoint const &A, LatticePoint const B, AuxDataType const &aux_data;
bool compare_abs_sc_product(LatticePoint const &A, LatticePoint const & B, ScalarProductReturnType target, AuxDataType const &aux_data); //compares whether the absolute value of the scalar product is at least target. This function might err.
bool compare_sc_product(LatticePoint const &A, LatticePoint const & B, ScalarProductReturnType target, AuxDataType const &aux_data); // Checks whether <A,B> > t


*/

template<class T> class DeclaresScalarProductReturnType
{
    private:
    template<class TT>
    static typename TT::ScalarProductReturnType foo(int);
    template<class ...>
    static void foo(...);
public:
    using value = std::is_void<decltype(foo<T>(0))>::value;
}

template<class Implementation>
class GeneralLatticePoint
{
    public:
    using LatticePointTag = true_type;
    using AuxDataType = IgnoreAnyArg;
    static_assert(DeclaresScalarProductReturnType<Implementation>::value, "Lattice Point class does not typedef its scalar product type");
    explicit constexpr GeneralLatticePoint()=default; //only ever called from its children
    GeneralLatticePoint(GeneralLatticePoint const &other)=delete; //This is just to match the implementation of a typical instantiation. Only the default constructor is ever used.
    GeneralLatticePoint(GeneralLatticePoint &&other)=default;
    GeneralLatticePoint& operator=(GeneralLatticePoint const & other) = delete;
    GeneralLatticePoint& operator=(GeneralLatticePoint && other) = default;
    ~GeneralLatticePoint()=default;
    typename Implementation::ScalarProductReturnType
    //The calling syntax for the derived object (i.e. Implementation) is supposed to be Implementation(DIM dim, AUX aux={}),
    //where DIM is convertible from int
    //and AUX is convertible from Implementation::AuxDataType
}



template<class ET,int nfixed>

    explicit MyLatticePoint(MatrixRow<ET> const & row, Dimension<nfixed> const & dim) {
        data = (row.get_underlying_row()).get();
        update_norm2( dim);
    };

    void update_norm2(Dimension<nfixed> const & dim)
    {
        this->norm2 = compute_sc_product(*this, *this, dim);
    }

    ET get_norm2() {return this->norm2;}

    //friend std::ostream & operator<< <ET, nfixed> (std::ostream &os, MyLatticePoint<ET,nfixed> const &A);


public:
    std::vector<ET> data;
    ET norm2;
};


template <class ET,int nfixed> MyLatticePoint<ET, nfixed> add (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> MyLatticePoint<ET,nfixed> minus (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET, nfixed> const &B, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> MyLatticePoint<ET,nfixed> negate (MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> void scalar_mult (MyLatticePoint <ET,nfixed> &A, ET const & multiple, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> bool comapre_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B,  ET target);
template <class ET,int nfixed> bool comapre_abs_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B,  ET target);
template <class ET,int nfixed> ET compute_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata);
template <class ET,int nfixed> MyLatticePoint<ET, nfixed> make_copy (MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> auxdata);


#endif
