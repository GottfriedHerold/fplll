#ifndef MY_LATTICE_POINT_CLASS_CPP
#define MY_LATTICE_POINT_CLASS_CPP


#include "MyLatticePointClass.h"


template <class ET,int nfixed>
MyLatticePoint<ET, nfixed> add (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata)
{
    MyLatticePoint<ET, nfixed> sum = MyLatticePoint<ET, nfixed>(auxdata.dim, auxdata);
    
    //MyLatticePoint<ET, nfixed> sum = make_copy(A, auxdata);
    
    for (unsigned int i=0; i<auxdata; ++i)
    {
        //sum.data[i] = sum.data[i]+B.data[i];
        sum.data[i].add(A.data[i], B.data[i]); //IS CORRECT?
    }
    
    sum.update_norm2(auxdata);
 
    return sum;
    
}

template <class ET, int nfixed>
MyLatticePoint<ET,nfixed> sub (MyLatticePoint<ET,nfixed> const &A, MyLatticePoint<ET, nfixed> const &B, Dimension<nfixed> const & auxdata)
{
    MyLatticePoint<ET, nfixed> sum = MyLatticePoint<ET, nfixed>(auxdata.dim, auxdata);
    
    for (unsigned int i=0; i<auxdata; ++i)
    {
        sum.data[i].sub(A.data[i], B.data[i]);
    }
    
    sum.update_norm2(auxdata);
    
    return sum;
}

template <class ET,int nfixed> MyLatticePoint<ET,nfixed> negateP (MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata)
{
    ET zero;
    zero = 0;
    MyLatticePoint<ET, nfixed> neg = MyLatticePoint<ET, nfixed>(auxdata.dim, auxdata);
    for (unsigned int i=0; i<auxdata; ++i)
    {
        neg.data[i].sub(zero, A.data[i]);
    }
    return neg;
}

template <class ET,int nfixed>
MyLatticePoint<ET,nfixed> scalar_mult (MyLatticePoint<ET,nfixed> &A, ET const & multiple, Dimension<nfixed> const & auxdata)
{
    MyLatticePoint<ET, nfixed> res = MyLatticePoint<ET, nfixed>(auxdata.dim, auxdata);
    
    for (unsigned int i=0; i<auxdata; ++i)
    {
        res.data[i].mul(A.data[i], multiple);
    }
    
    res.update_norm2(auxdata);
    
    return res;

    
}

template <class ET,int nfixed>
bool compare_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B,  ET target)
{
    return false;
}


template <class ET,int nfixed>
ET compute_sc_product (MyLatticePoint<ET, nfixed> const &A, MyLatticePoint<ET,nfixed> const &B, Dimension<nfixed> const & auxdata)
{
    ET res;
    res = 0;
    if (auxdata>0)
    {
        res.mul(A.data[0],B.data[0]);
        for (unsigned int i=1; i<auxdata; i++)
        {
            res.addmul(A.data[i], B.data[i]);
        }
        
    }
    
    return res;
}


template <class ET,int nfixed>
MyLatticePoint<ET, nfixed> make_copy (MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata)
{
    MyLatticePoint<ET, nfixed> copy (A, auxdata);
    return copy;
}

/*
template <class ET,int nfixed> MyLatticePoint<ET, nfixed>
void print (std::ostream &os, MyLatticePoint<ET,nfixed> const &A, Dimension<nfixed> const & auxdata)
{
    
}
*/
#endif
