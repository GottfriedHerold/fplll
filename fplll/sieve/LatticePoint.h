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

 */


template<class ZT>
class LatticePoint {

public:

    /* vector */
    NumVect<Z_NR<ZT> > v;

    /* square L2 norm of the vector */
    Z_NR<ZT> norm;

    LatticePoint<ZT> (int n)
    {
       norm = 0;
       v.resize(n);
        for (int i=0; i<n; i++)
        {
           v[i] = 0;
        }

    }

};

template <class ZT> inline LatticePoint<ZT>* new_latticepoint(int n)
{
    LatticePoint<ZT> *p = new LatticePoint<ZT> (10);
    p->norm = 0;
    p->v.resize(n);
    for (int i = 0; i < n; ++i)
    {
        p->v[i] = 0;
    }
    return p;
}


template <class ZT>
void printLatticePoint (LatticePoint<ZT> *p)
{
    cout << p->v << endl;
}



//template<class ZT>
//class LatticePoint : protected NumVect<ZT>
//{
//public:
//    /* square L2 norm of the vector */
//        Z_NR<ZT> norm;
//
//}
#endif
