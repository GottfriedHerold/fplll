#ifndef GAUSS_SIEVE_UTILITY_H
#define GAUSS_SIEVE_UTILITY_H

namespace GaussSieve //helper functions
{
    bool string_consume(istream &is, std::string const & str, bool elim_ws= true, bool verbose=true); //helper function for dumping/reading
    Z_NR<mpz_t> compute_mink_bound(ZZ_mat<mpz_t> const & basis);
    template<class ET>
    bool check_perform_2red (LatticePoint<ET> &p1, const LatticePoint<ET> &p2); //2-reduces p1 with the help of p2.
                                                                       //p1 is overwritten, whereas p2 is const. Returns true if p1 actually changed.
    template<class ET>
    bool check2red_new (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET &scalar); //only checks whether 2reduction is possible

    template<class ET>
    LatticePoint<ET> perform2red (const LatticePoint<ET> &p1, const LatticePoint<ET> &p2, ET const & scalar); //replaces p1 by p1 - scalar * p2

    template<class ET>
    bool check3red(const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, float px1, float px2, float x1x2, int & sgn1, int & sgn2);

    template<class ET>
    bool check3red_signs(const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, int px1, int px2, int x1x2, int &sgn1, int &sgn2);

    template<class ET>
    LatticePoint<ET> perform3red (const LatticePoint<ET> &p, const LatticePoint<ET> &x1, const LatticePoint<ET> &x2, const int & sgn1, const int &sgn2);
}

#include "Utility.cpp"
#endif // GAUSS_SIEVE_UTILITY_H
