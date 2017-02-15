//
//  LatticePoint.h
//
#ifndef LATTICE_VECTOR_CLASS_H
#define LATTICE_VECTOR_CLASS_H

#include "sieve_common.h" //needed (at least for convertions from MatrixRow (the header has to be revised);

/**
 * Class for list element

    data types:

    class Z_NR: stores integers; defined in type_nr_Z.h
    class NumVect: uses STL's container of type <vector>; defined in type_numvec.h


    members:
    v: instance of NumVect
    norm: l^2-norm of v

**/

template <class ET> class LatticePoint;
template <class ET> class IsLongerVector_class; //class wrapper to compare vectors by length. Needs to be wrapped in a class to work seamlessly with some STL containers.

template <class ET>
void sc_product (ET &result, const LatticePoint<ET> &p1, const LatticePoint<ET> &p2);



template<class ET> //ET: entries of the individual vectors. Required to be copy-constructible. Use E = Z_NR<mpz_t> rather than E=mpz_t.
class LatticePoint : public NumVect<ET>
{
       using NV = NumVect<ET>;

    /* square L2 norm of the vector */
public:
        ET norm2;

public:
    LatticePoint(){}
    LatticePoint(const LatticePoint &Point) = default; // : NumVect<ET>::data(Point.data), norm2(Point.norm2) {}
    LatticePoint(LatticePoint &&Point) = default ; // : NumVect<ET> (), norm2(0) {swap(Point.NV::data), swap(Point.norm2);}
    LatticePoint(int n) : NumVect<ET>(n), norm2() //creates all-zero vector of length n
        {
           this->data.resize(n);
           this->fill(0);
	       norm2 = 0; //why? LatticePoint(int n) : NumVect<ET>(n), norm(0) errs
        }

    // for debugging
    LatticePoint(int n, long fillwith) : NumVect<ET>(n),norm2()
    {
        this->data.resize(n);
        this->fill(fillwith);
        sc_product(this->norm2, *this, *this);
        //norm2 = fillwith*fillwith*n; gives type problems.
    }
    //TODO: fillwith should be ET, not long.

    LatticePoint(NumVect<ET> vector_) : NumVect<ET>(vector_)
    {
    sc_product(this->norm2, *this, *this);
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
    inline ET get_norm2() const {return norm2;}
    inline void get_norm2 (ET norm_to_return) const {norm_to_return = norm2;}

    inline void setNorm2 (ET norm) {this->norm2 = norm;}

    void printLatticePoint()
    {
        cout << * (static_cast<NumVect<ET>*>(this)) << " of norm: " << this->norm2 << endl;
    }
//    void subtr (LatticePoint<ET> const A)
//    {
//        this->this.sub(A, A.size());
//    }
};

template<class ET> class IsLongerVector_class
{
    public: bool operator() (LatticePoint<ET> const &A, LatticePoint<ET> const & B)
    {
     return (A.get_norm2() > B.get_norm2() );
    }
};


template <class ET>
LatticePoint<ET> operator+ (LatticePoint<ET> const &A, LatticePoint<ET> const &B)
{

	LatticePoint<ET> C(A);
	//length-check is done in by add in numvect
	C.add(B, A.size());
	ET norm;
	sc_product (norm, C, C);
	C.setNorm2(norm);
	return C;

}

template <class ET>
LatticePoint<ET> operator- (LatticePoint<ET> const &A, LatticePoint<ET> const &B)
{

	LatticePoint<ET> C(A);
	C.sub(B, A.size());
	ET norm;
	sc_product (norm, C, C);
	C.setNorm2(norm);
	return C;

}

template <class ET>
void scalar_mult (LatticePoint<ET> &A, ET multiple)
{
    A.mul(A, A.size(), multiple);
}


//Simple dot_product
template <class ET>
void sc_product (ET &result, const LatticePoint<ET> &p1, const LatticePoint<ET> &p2)
{
   dot_product(result, p1, p2);
}

//Convert MatrixRow to LatticePoint

template <class ET>
LatticePoint<ET> conv_matrixrow_to_lattice_point (MatrixRow<ET> const &row)
{
	LatticePoint<ET> res(row.get_underlying_row());
	//NumVect<ET> tmp(row.get_underlying_row());
	return res;
}

// Convert sample() result NumVect to LatticePoint

template<class ET>
inline LatticePoint<ET> conv_sample_to_lattice_point (NumVect<ET> const &sample)
{
    LatticePoint<ET> p (sample);
    return p;
}

#endif
