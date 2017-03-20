#ifndef ELLIPTIC_SAMPLER_H
#define ELLIPTIC_SAMPLER_H

#include "Sampler.h"

#include <random>

template<class ET, bool MT, class Sseq = std::seed_seq>
class EllipticSampler;



template<class ET,bool MT,  class Sseq>
class EllipticSampler: public Sampler<ET,MT, Sseq>
{
    public:
    EllipticSampler(double const s=2.0, double const cutoff = 2.0) : Sampler<ET,MT,Sseq>(),s2pi(s*s*GaussSieve::pi),maxdeviations(s*cutoff) {};
    virtual void init(Sieve<ET,MT> * const sieve, Sseq seed) override;
    virtual SamplerType  sampler_type() const override                          {return SamplerType::elliptic_sampler;};
    virtual ~EllipticSampler();
    virtual LatticePoint<ET> sample(int thread=0) override;
    private:
    Sieve<ET,MT> * sieveptr; //pointer to parent sieve.
    ZZ_mat<typename ET::underlying_data_type> current_basis;
    Matrix<FP_NR<double > > mu;
    double s2pi; //stores standard dev. for each dimension, already squared and multiplied by pi.
    double maxdeviations; //stores s*cutoff for each dimension.
    unsigned int rank;

};

template<class ET,bool MT, class Sseq>
void EllipticSampler<ET,MT,Sseq>::init(Sieve<ET,MT> * const sieve, Sseq seed)
{
    sieveptr = sieve;
    current_basis = sieve->get_original_basis();
    rank = sieve->get_lattice_rank();
    Matrix<ET> u, u_inv,g; //intentionally uninitialized.
    MatGSO<ET, FP_NR<double> > GSO(current_basis, u, u_inv, MatGSOFlags::GSO_INT_GRAM);
    GSO.update_gso(); //todo: raise exception in case of error.
    mu = GSO.get_mu_matrix();
//    g  = pGSO.get_g_matrix();
//    F maxbistar2 = pGSO.get_max_bstar();
//    F tmp;
//    variances.resize(rank); variances.shrink_to_fit();
//    maxdeviations.resize(rank); maxdeviations.shrink_to_fit();
//

}
template<class ET,bool MT, class Sseq>
EllipticSampler<ET,MT,Sseq>::~EllipticSampler()
{

}

template<class ET,bool MT, class Sseq>
LatticePoint<ET> EllipticSampler<ET,MT,Sseq>::sample(int thread)
{

}



#endif
