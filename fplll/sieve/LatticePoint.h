//
//  LatticePoint.h
//
#ifndef LATTICE_VECTOR_CLASS_H
#define LATTICE_VECTOR_CLASS_H

#include "sieve_common.h"

/**
 * Class for list element

    data types:

    class Z_NR: stores integers; defined in type_nr_Z.h
    class NumVect: uses STL's container of type <vector>; defined in type_numvec.h


    members:
    v: instance of NumVect
    norm: l^2-norm of v

**/

template<class ZT>
class LatticePoint : public NumVect<ZT>
{
    /* square L2 norm of the vector */
        ZT norm2;

   using NV = NumVect<ZT>;

public:
    LatticePoint(){}
    LatticePoint(const LatticePoint &Point) = default; // : NumVect<ZT>::data(Point.data), norm2(Point.norm2) {}
    LatticePoint(LatticePoint &&Point) = default ; // : NumVect<ZT> (), norm2(0) {swap(Point.NV::data), swap(Point.norm2);}
    LatticePoint(int n) : NumVect<ZT>(n), norm2(0) //creates all-zero vector of length n
        {
           this->data.resize(n);
           this->fill(0);
	   //norm2 = 0; //why? LatticePoint(int n) : NumVect<ZT>(n), norm(0) errs
        }

    // for debugging
    LatticePoint(int n, long fillwith) : NumVect<ZT>(n),norm2()
    {
        this->data.resize(n);
        this->fill(fillwith);
	
    }

    LatticePoint& operator=(LatticePoint that)
        {
            swap(that);
            return *this;
        }

    ~LatticePoint() {}

    void swap(LatticePoint &Point)
    {
        this->data.swap(Point.data);
        std::swap(norm2, Point.norm2);
    }

    inline NV& getVector() const {return this->data.get();}
    inline ZT getNorm() const {return norm2;}

    inline void setNorm2 (ZT norm) {this->norm2 = norm;}

    void printLatticePoint()
{
    //using NV = NumVect<ZT>;
    cout << *this << " of norm: " << this->norm2 << endl;
    //cout << this->getNorm() << endl;
}


};

template<class ZT>
class IsLongerVector_class //should be moved to LatticePoint.h. Make sure getNorm is declared const.
{
    public:
    bool operator() (LatticePoint<ZT> const &A, LatticePoint<ZT> const & B)
    {
     return (A.getNorm() > B.getNorm() );
    }
};

template <class ZT>
LatticePoint<ZT> operator+ (LatticePoint<ZT> const &A, LatticePoint<ZT> const &B)
{	

	LatticePoint<ZT> C(A);
	//length-check is done in by add in numvect
	C.add(B, A.size());
	ZT norm;
	sc_product (norm, C, C);
	C.setNorm2(norm);
	return C;

}

//unary minus
template <class ZT>
LatticePoint<ZT> operator- (LatticePoint<ZT> const &A, LatticePoint<ZT> const &B)
{	

	LatticePoint<ZT> C(A);
	C.sub(B, A.size());
	ZT norm;
	sc_product (norm, C, C);
	C.setNorm2(norm);
	return C;

}


//Simple dot_product
template <class ZT>
void sc_product (ZT &result, const LatticePoint<ZT> &p1, const LatticePoint<ZT> &p2)
{
    //using NV = NumVect<ZT>;
   dot_product(result, p1, p2);
}

#endif
