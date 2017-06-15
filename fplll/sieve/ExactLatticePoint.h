#ifndef EXACT_LATTICE_POINT_H
#define EXACT_LATTICE_POINT_H

/*
    ExactLatticePoint stores an (exact) lattice point, together with its squared length.
    We also have some elementary arithmetic (addition, substraction etc. ) defined on them.
*/


template<class ET,int nfixed> //ET: entries of the individual vectors. Required to be copy-constructible. Use E = Z_NR<mpz_t> rather than E=mpz_t.
class ExactLatticePoint : public NumVect<ET>
{
       using NV = NumVect<ET>;

    /* square L2 norm of the vector */
public:
        ET norm2;
public:
    ExactLatticePoint()=default;
    ExactLatticePoint(const ExactLatticePoint &Point) = default;
    ExactLatticePoint(ExactLatticePoint &&Point) = default ;
    ExactLatticePoint(int n); //creates all-zero vector of length n
    ExactLatticePoint(int n, long fillwith);    //TODO: fillwith should be ET, not long.
    ExactLatticePoint(NumVect<ET> vector_) : NumVect<ET>(vector_) {normalize();} //Creates an exact lattice point from a NumVect
    void normalize(); //sets norm2 to the correct value. Use after using operations from the NumVect<ET> base class.
    ExactLatticePoint& operator=(ExactLatticePoint const &that) =default;
    ExactLatticePoint& operator=(ExactLatticePoint && that) =default;
    ~ExactLatticePoint() {}

    friend void swap(ExactLatticePoint &A, ExactLatticePoint &B)
    {
        std::swap(A.data, B.data);
        std::swap(A.norm2,B.norm2);
    }
    friend ostream & operator<<(ostream &os, ExactLatticePoint const & exact_point) {exact_point.print_lattice_point(os); return os;} //printing
    //friend istream & operator>>(istream &is, ExactLatticePoint & exact_point); //reading (may also be used by constructor from istream)

    NV const & access_vector() const {return *this;} //Probably not needed anyway.
    ET const & access_norm2() const {return norm2;} //Note: Does not copy
    unsigned int get_dim() const {return NumVect<ET>::size();} //returns dimension
    //inline void get_norm2 (ET &norm_to_return) const {norm_to_return = norm2;}
    //inline void setNorm2 (ET norm) {this->norm2 = norm;} //should not be required, actually.

    void print_lattice_point(ostream & os = cout) const
    {
        os << * (static_cast<NumVect<ET> const *>(this)) << " of norm: " << this->norm2 << endl;
    }
};

#include "ExactLatticePoint.cpp"
#endif
