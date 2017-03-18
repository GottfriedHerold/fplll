#ifndef SAMPLER_H
#define SAMPLER_H

//forward declarations
#include <random>
#include <iostream>
#include <cfenv>

template<class ET, bool MT, class Sseq = std::seed_seq>
class Sampler;
template<class ET,bool MT> ostream & operator<<(ostream &os, Sampler<ET,MT>* const samplerptr); //printing
template<class ET,bool MT> istream & operator>>(istream &is, Sampler<ET,MT>* const samplerptr); //reading (may also be used by constructor from istream)
enum class SamplerType
{
    user_defined = 0,
    klein_old    = 1
};
template<class Engine, bool MT, class Sseq = std::seed_seq> //make separate class to allow specialisation for MT.
class MTPRNG;                     //wrapper around (a vector of) random number engines of type Engine
template<class ET, bool MT, class Sseq = std::seed_seq>
class KleinSamplerNew;
template<class ET, bool MT, class F, class Sseq = std::seed_seq>
class KleinSamplerOld;

namespace GaussSieve
{
    template<class Z, class Engine>
    Z sample_z_gaussian(double s, double center, Engine & engine, double cutoff);
    long double constexpr pi_long = std::atan2(0, -1); // pi = 3.14...
    double constexpr pi_double = std::atan2(0,-1);
    long double constexpr pi=std::atan2(0,-1);
    //samples from a discrete Gaussian distribution with parameter s and center c. We cutoff the Gaussian at s*cutoff.
    //i.e. the distribution is discrete on Z with output probability for x being proportional to exp(- pi(x-c)^2/s^2). Note the scaling by pi in the exponent.
    //For reasons of numerical stability, center should not be very large in absolute value (it is possible to reduce to |center|<1 anyway), s.t.
    //center +/- cutoff * s does not overflow.
    //Z must be among one of short, int, long, long long.
    //We do NOT support mpz_t here! Note that the output takes the role of coefficients wrt a given basis.
    //We only support double. For sieving algorithms, there is no really good reason to support higher precision.
    //Note: if one wants to have higher precision, one also needs to adjust the PRNGs to actually output high precision.

    //inline double gaussian_fn(double s, double x);                  //computes the (scaled) Gaussian density exp(-pi x^2 / s^2)
    //inline double gaussian_fn_adj(double s, double x, double adj);  //computes the (scaled) Gaussian density exp(-pi (x^2 - adj)/s^2 )
    //Note: the version with adj is proportional to the one without, so it represents the same distribution. It's intended use is as follows:
    //If x = z+c with z integral, set adj to be min( (z+c)^2). This way, the density is scaled up (and reaches 1 at the most likely value).
    //If s<<dist(c,\mathbb{Z})<=1/2, the gaussian_fn variant always outputs extremely small values, leading to bad performance in rejection sampling.
    //The adjusted version avoids this problem.
}

//includes
#include "SieveGauss.h"


//The class MTPRNG is just a wrapper around a PRNG engine to facilitate switching to multi-threaded. Due to the fact that we support multi-threading, MTPRNG<Engine,true,.> is a wrapper around
//a vector of Engines, whereas, MTPRNG<Engine,false,.> is a wrapper around Engine.
//reseed seeds *all* engines. Use rnd(thread-id) to obtain the underlying Engine.

template<class Engine, class Sseq>
class MTPRNG<Engine,true, Sseq>
{
    public:
    MTPRNG(unsigned int const _num_threads) : engines(_num_threads), num_threads(_num_threads) {};
    void reseed(Sseq & seq);
    Engine & rnd(int const thread)            {return engines[thread];};
    private:
    std::vector<Engine> engines;
    unsigned int const num_threads;
};

template<class Engine, class Sseq>
void MTPRNG<Engine,true,Sseq>::reseed(Sseq & seq)
{
    std::vector<typename Sseq::result_type> seeds(num_threads);
    seq.generate(seeds.begin(), seeds.end());
    for(int i=0;i<num_threads;++i)
    {
        engines[i].seed(seeds[i]);
    }
}

template<class Engine, class Sseq> //just wrapper around Engine
class MTPRNG<Engine, false, Sseq>
{
    public:
    MTPRNG(int const num_threads = 1) : engine() {};
    void reseed(Sseq & seq)                 {engine.seed(seq);};
    Engine & rnd(int const thread = 1)      {return engine;};
    private:
    Engine engine;
};      //End of MTPRNG

//generic Sampler. All other sampler are derived from it.

template<class ET,bool MT, class Sseq> //Sseq is supposed to satisfy the C++ concept "SeedSequence". The standard library has std::seed_seq as a canonical example.
class Sampler
//Note :    In multi-threaded environment, we only have 1 sampler object
//          caller_thread is set to -1 if we make a single-threaded call
{
    public:
    friend ostream & operator<< <ET,MT>(ostream &os, Sampler<ET,MT,Sseq>* const samplerptr);
    friend istream & operator>> <ET,MT>(istream &is, Sampler<ET,MT,Sseq>* const samplerptr);

    virtual void init(Sieve<ET,MT> * const sieve, Sseq seed) {}                           //called before any points are sampled;
    virtual ~Sampler()=0; //needs to be virtual
    virtual SamplerType  sampler_type() const {return SamplerType::user_defined;};    //run-time type information.
                                                                    //This may be used to determine how to interpret a dump file.
                                                                    //defaults to user-defined.
                                                                    //Other values mean that the GaussSieve dumping routine is aware of the type, simplifying the syntax for dumping / reading.
    //TODO : Write some explanation how to do that.
    virtual LatticePoint<ET> sample(int thread=0)=0; //thread is the index of the calling thread (we need to keep separate PRNGs for each thread)
    //TODO : Allow sampling in subspaces, updating basis.

    private:
    virtual ostream & dump_to_stream(ostream &os)  {return os;};    //dummy implementation of << operator.
    virtual istream & read_from_stream(istream &is){return is;};    //dummy implementation of >> operator.

};
template <class ET,bool MT, class Sseq>
Sampler<ET,MT, Sseq>::~Sampler() {} //actually needed, even though destructor is pure virtual as the base class destructor is eventually called implicitly.

template<class ET,bool MT> ostream & operator<<(ostream &os,Sampler<ET,MT>* const samplerptr){return samplerptr->dump_to_stream(os);};
template<class ET,bool MT> istream & operator>>(istream &is,Sampler<ET,MT>* const samplerptr){return samplerptr->read_from_stream(is);};


//GPV / Klein sampler, adapted from old sampler_basic.*, but using our own conventions.

template<class ET,bool MT, class F,class Sseq>
class KleinSamplerOld : public Sampler<ET,MT, Sseq>
{
    public:
    virtual void init(Sieve<ET,MT> * const sieve, Sseq seed) override;
    virtual SamplerType  sampler_type() const override                          {return SamplerType::klein_old;};
    virtual ~KleinSamplerOld();
    virtual LatticePoint<ET> sample(int thread=0) override;
    private:
    Sieve<ET,MT> * sieveptr; //pointer to parent sieve.
    ZZ_mat<typename ET::underlying_data_type> current_basis;
};

template<class ET,bool MT, class F, class Sseq>
void KleinSamplerOld<ET,MT,F,Sseq>::init(Sieve<ET,MT> * const sieve, Sseq seed)
{

}
template<class ET,bool MT, class F, class Sseq>
KleinSamplerOld<ET,MT,F,Sseq>::~KleinSamplerOld()
{

}

template<class ET,bool MT, class F, class Sseq>
LatticePoint<ET> KleinSamplerOld<ET,MT,F,Sseq>::sample(int thread)
{

}

/**
 * sampling Z by rejection sampling
 */

template<class Z, class Engine>
Z GaussSieve::sample_z_gaussian(double s, double center, Engine & engine, double cutoff)
{
//Note : The following allows to access / modify floating point exceptions and modes.
//#pragma STDC FENV_ACCESS on
//This is too compiler/implementation-specific and does not work most of the time...

    static_assert(is_integral<Z>::value,"Return type for sample_z_gaussian must be POD integral type.");
    Z maxdev = static_cast<Z>(std::ceil(s * cutoff)); //maximum deviation of the Gaussian from the center. Note that maxdev may be 1.
    std::uniform_int_distribution<Z> uniform_in_range (std::floor(center-maxdev),std::ceil(center+maxdev));
    std::uniform_real_distribution<double> rejection_test(); //defaults to value from [0,1), used in rejection sampling.
    Z closest_int = std::round(center); //closest int to center, i.e. most likely value.
    double adj = -(center-closest_int)*(center-closest_int); //negative squared distance to most likely value. Used to scale up the Gaussian weight function s.t. it is 1 at the most likely value.

    s = s*s/pi; //overwriting s.
    //std::fenv_t env;
    //feholdexcept( &env); //This disables all floating-point exceptions.

//use rejection sampling
    while(true)
    {
        Z result = uniform_in_range(engine); //sample uniform result.
        double dist = result - center;
    //compute Gaussian weight. std::fma(dist,dist,adj) computes dist^2 + adj = (result-center)^2  - MIN{(result-center)^2 | result integral}.
    //s was overwritten to be s^2/pi.

    //Note that the argument of the exp-function might be a tiny positive value due to numeric error
    //(Even if result==closest_int, adj = ROUND((closest_int-center)^2), the computation of std::fma(dist,dist,adj) does not round the intermediate dist^2, leading to a non-zero argument)
    //In particular, it is conceivable that floating point underruns occur in the std::fma - call.
    //Furthermore, if cutoff is large or if s<<1 (in this case, the issue is the rounding when we determined the range), the argument to exp can be extremely small, leading to further potential underruns.
    //We do not care about this for now...

        if( rejection_test(engine) <  std::exp(-std::fma(dist,dist,adj)/s))
        {
            //std::feclearexcept(FE_UNDERFLOW);
            //std::feupdateenv(&env);
            return result;
        }

}
}
//Old code:
//  F min, max, st, range, tmp, tmp1;
//  double r, e;
//
//  /* (c \pm s*t) for t \approx logn */
//  st = s;
//  st.mul(st, t, GMP_RNDN);
//  min.sub(c, st, GMP_RNDN);
//  max.add(c, st, GMP_RNDN);
//  min.rnd(min);
//  max.rnd(max);
//  range.sub(max, min, GMP_RNDN);
//
//  Z_NR<ZT> x;
//  while (1)
//  {
//    r = double(rand()) / RAND_MAX;
//    tmp.mul_d(range, r, GMP_RNDN);
//    tmp.rnd(tmp);
//    tmp.add(tmp, min, GMP_RNDN);
//    x.set_f(tmp);
//    tmp1.sub(tmp, c, GMP_RNDN);
//    tmp1.mul(tmp1, tmp1, GMP_RNDN);
//    tmp1.mul_d(tmp1, -M_PI, GMP_RNDN);
//    tmp.mul(s, s, GMP_RNDN);
//    tmp1.div(tmp1, tmp, GMP_RNDN);
//    e = tmp1.get_d(GMP_RNDN);
//    r = exp(e);
//    if ((double(rand()) / RAND_MAX) <= r)
//      return x;
//  }


#endif
