#ifndef LATTICE_POINT_CLASS_CPP
#define LATTICE_POINT_CLASS_CPP


#include "LatticePointClass.h"


template <class ET,int nfixed>
MyLatticePoint<ET> add (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET,nfixed> const &B)
{
    
    MyLatticePoint<ET,nfixed> C(A);
    return C;
    
}

template <class ET, int nfixed>
MyLatticePoint<ET,nfixed> minus (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET> const &B)
{
    
    MyLatticePoint<ET,nfixed> C(A);
    return C;
}

template <class ET,int nfixed> MyLatticePoint<ET,nfixed> negate (MyLatticePoint<ET,nfixed> const &A)
{
    MyLatticePoint<ET,nfixed> neg(A);
    return neg;
}

template <class ET,int nfixed>
void scalar_mult (MyLatticePoint<ET,nfixed> &A, ET const & multiple)
{
    
}

template <class ET,int nfixed>
bool comapre_sc_product (MyLatticePoint<ET, fixed> const &A, MyLatticePoint<ET,nfixed> const &B,  ET target)
{
    return false;
}

template <class ET,int nfixed>
ET compute_sc_product (MyLatticePoint<ET, fixed> const &A, MyLatticePoint<ET,nfixed> const &B)
{
    
}

template <class ET,int nfixed>
MyLatticePoint<ET, nfixed> make_copy (MyLatticePoint<ET,nfixed> const &A, AuxDataType auxdata)
{
    
}

#endif
