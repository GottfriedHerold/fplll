#include "ShiSampler.h"
#include "Sampler.h"
#include "LatticePointsNew.h"

template<class ET,bool MT, class Engine, class Sseq> void ShiSampler<ET,MT,Engine,Sseq,-1>::custom_init()
{
    current_basis = sieveptr->get_original_basis();
    dim = sieveptr->get_ambient_dimension();
    rank = sieveptr->get_lattice_rank();
    Matrix<ET> u, u_inv,g; //intentionally uninitialized.
    MatGSO<ET, FP_NR<double> > GSO(current_basis, u, u_inv, MatGSOFlags::GSO_INT_GRAM);
    GSO.update_gso(); //todo: raise exception in case of error.
    mu = GSO.get_mu_matrix();

    s2pi.resize(rank);
    maxdeviations.resize(rank);
    g  = GSO.get_g_matrix();

    FP_NR<double> maxbistar2 = GSO.get_max_bstar();
    FP_NR<double> tmp;
    FP_NR<double> tmp2;
    for (unsigned int i = 0; i < rank; ++i)
    {
        tmp.set_z(g(i, i));
        tmp2.div(maxbistar2,tmp); //s'_i^2 = max GS length^2 / lenght^2 of basis vector
        s2pi[i] = tmp2.get_d() / GaussSieve::pi;
        tmp2.sqrt(tmp2);
        maxdeviations[i] = tmp2.get_d() * cutoff;
    }
}

template<class ET,bool MT, class Engine, class Sseq> ShiSampler<ET,MT,Engine, Sseq>::~ShiSampler()
{

}

template<class ET,bool MT, class Engine, class Sseq> typename  Sampler<ET,MT,Engine,Sseq,-1>::SampleReturnType ShiSampler<ET,MT,Engine, Sseq>::sample(int thread)
{
    assert(sieveptr!=nullptr);
    ExactLatticePoint<ET,-1> *vec = new ExactLatticePoint<ET,-1>(dim);
    vec->NumVect<ET>::fill(0); //current vector built up so far. //Note: We treat vec as a NumVect until the end, because we don't want to normalize intermediate results.
    vector<double> shifts(rank, 0.0); //shift, expressed in coordinates wrt the Gram-Schmidt basis.
    for(int j=rank-1; j>=0; --j)
    {
        long const newcoeff = GaussSieve::sample_z_gaussian_VMD<long,Engine>(s2pi[j],shifts[j],engine.rnd(),maxdeviations[j]); //coefficient of b_j in vec.
        //vec+= current_basis[j].get_underlying_row(); //build up vector
        vec->NumVect<ET>::addmul_si(current_basis[j].get_underlying_row(), newcoeff);
        for(int i=0;i<j;++i) //adjust shifts
        {
            shifts[i]-=newcoeff* (mu[j][i].get_d() );
        }
    }
    vec->normalize();

    return static_cast<CompressedPoint<ET,MT,-1> > (vec); //Note: This makes copies, which are unneccessary...
}
