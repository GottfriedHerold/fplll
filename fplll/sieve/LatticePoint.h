//
//  LatticePoint.h
//
//

#ifndef _lattice_vector_class_h
#define _lattice_vector_class_h

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
class LatticePoint : protected NumVect<ZT>
{
    /* square L2 norm of the vector */
        Z_NR<ZT> norm2;

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
    LatticePoint(int n, long fillwith) : NumVect<ZT>(n),norm2(0)
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

    inline NV& getVector() {return this->data.get();}
    inline Z_NR<ZT> getNorm() {return norm2;}

};

template <class ZT>
void printLatticePoint (LatticePoint<ZT> &p)
{
    //using NV = NumVect<ZT>;
    cout << p.getNorm() << endl;
}


//Simple dot_product
template <class ZT>
void dot_product (ZT &result, const LatticePoint<ZT> &p1, const LatticePoint<ZT> &p2)
{
    //using NV = NumVect<ZT>;
    dot_product(result, p1.data, p2.data);
}

#endif
