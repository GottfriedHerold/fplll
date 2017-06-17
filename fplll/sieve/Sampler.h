#ifndef SAMPLER_H
#define SAMPLER_H

//forward declarations
#include <random>
#include <iostream>
#include <cfenv>
#include <type_traits>

//forward declarations

template<class ET,bool MT, class Engine, class Sseq> class Sampler;
template<class ET,bool MT, class Engine, class Sseq> ostream & operator<<(ostream &os, Sampler<ET,MT, Engine, Sseq>* const samplerptr); //printing
template<class ET,bool MT, class Engine, class Sseq> istream & operator>>(istream &is, Sampler<ET,MT, Engine, Sseq>* const samplerptr); //reading (may also be used by constructor from istream)
enum class SamplerType
{
    user_defined = 0,
    elliptic_sampler= 1,
    shi_sampler = 2,
    gauss_sampler =3
};

template<class Engine, bool MT, class Sseq> class MTPRNG;       //wrapper around (a vector of) random number engines of type Engine
                                                                //This is used to unify the single and multi-threaded case
namespace GaussSieve
{
//    #if __GNUG__ //GCC has constexpr variants of trigonometric functions. Unfortunately, other compilers also define __GNUC__ and so on and detecting this at compile time is a pain.
//    long double constexpr pi_long = std::atan2(0, -1); // pi = 3.14...
//    double constexpr pi_double = std::atan2(0,-1);
//    long double constexpr pi=std::atan2(0,-1);
//    #else
    long double constexpr pi_long = 3.14159265358979323846264338327950288419716939937510L;
    double constexpr pi_double = 3.14159265358979323846264338327950288419716939937510;
    long double constexpr pi = 3.14159265358979323846264338327950288419716939937510L;
//    #endif // __GNUC__

    template<class Z, class Engine>
    Z sample_z_gaussian(double s, double const center, Engine & engine, double const cutoff);
    template<class Z, class Engine>
    Z sample_z_gaussian_VMD(double const s2pi, double const center, Engine & engine, double const maxdeviation);
    //samples from a discrete Gaussian distribution with parameter s and center c. We cutoff the Gaussian at s*cutoff.
    //i.e. the distribution is discrete on Z with output probability for x being proportional to exp(- pi(x-c)^2/s^2). Note the scaling by pi in the exponent.
    //For reasons of numerical stability, center should not be very large in absolute value (it is possible to reduce to |center|<1 anyway), s.t.
    //center +/- cutoff * s does not overflow.
    //Z must be among one of short, int, long, long long. center needs to be representable as an (exact) double.
    //We do NOT support mpz_t here! Note that the output takes the role of coefficients wrt a given basis.
    //We only support double. For sieving algorithms, there is no really good reason to support higher precision.
    //Note: if one wants to have higher precision, one also needs to adjust the PRNGs to actually output high precision.

    //The variant sample_z_gaussian_VMD takes s2pi = s^2 * pi and cutoff * s as parameters.


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
//Thread-safety: Init and reseed are not thread-safe. Concurrent calls to rnd are fine.
//You may concurrently call rnd and use init to increase the number of threads.

//Randomness path:
//The global (master) seed that is input to the constructor resp. reseed is used to create 20x32 bit per-thread-seeds for each actual Engine.
//The output from theses engine(s) is then accessed using rnd(thread-number).
//This is done even in the single-threaded case to ensure consistency.
//Note that for obtaining the per-thread seeds from the master seeds, we use a fixed Mersenne twister engine and not the template param.

template<class Engine, class Sseq>  class MTPRNG<Engine,true, Sseq>
{
    public:
    MTPRNG(Sseq & _seq = {}) : seeder(_seq), engines(0), num_threads(0)       {}; //constructs an uninitialized MTPRNG
    void reseed(Sseq & _seq);
    void init(unsigned int const _num_threads); //will make sure at least _num_threads engines are actually running, starting new ones as desired.
                                                //Will never reseed/restart already running engines.
                                                //Reducing the number of threads and increasing it back saves the random state (unless we reseed).
    Engine & rnd(unsigned int const thread)                                 {return engines[thread];};
    private:
    std::mt19937 seeder; //seeded with initial seq and consecutively used to seed the children PRNGs.
    std::vector<Engine> engines;
    unsigned int num_threads; //number of initialized engines. May differ from size of the vector. In particular, num_threads = 0 means uninitialized.
    Sseq seq;
    static unsigned int constexpr seed_length = 20; //number of 32bit values to use as seeds for the underlying engine(s).
                                                    //Technically, we could use state_size if the engine provides it, but not all default engines do.
};


template<class Engine, class Sseq>  class MTPRNG<Engine, false, Sseq>           //just wrapper around Engine
{
    public:
    MTPRNG(Sseq & _seq ={}) : engine()                      {reseed(_seq);};
    void reseed(Sseq & _seq);
    void init(unsigned int const = 1)                       {} //does nothing.
    Engine & rnd(unsigned int const)                        {return engine;};   //Argument is number of thread. It is ignored.
    Engine & rnd()                                          {return engine;};   //Version without thread-id
    private:
    Engine engine;
    static unsigned int constexpr seed_length = 20; //number of 32bit values to use as (per-thread) seed for the underlying engine.
                                                    //Technically, we could use state_size if the engine provides it, but not all default engines do.
};      //End of MTPRNG

template<class Engine, class Sseq> void MTPRNG<Engine,true,Sseq>::reseed(Sseq & _seq)
{
    seeder.seed(_seq);
    unsigned int old_threads=num_threads;
    num_threads=0;
    init(old_threads); //will restart all engines, because num_threads = 0;
};

template<class Engine, class Sseq> void MTPRNG<Engine,true,Sseq>::init(unsigned int const _num_threads)
{
    if(_num_threads<=num_threads) //no need to initalize.
    {
        return;
    }
    engines.resize(_num_threads);
    engines.shrink_to_fit();
    uint32_t per_engine_seed[seed_length];
    //else initialize remaining threads
    for(int i=num_threads;i<_num_threads;++i)
    {

        for(int j=0;j<seed_length;++j)
        {
            per_engine_seed[j] = seeder();
        }
        std::seed_seq per_engine_see_seq(per_engine_seed, per_engine_seed+seed_length);
        engines[i].seed(per_engine_see_seq);
    }
    num_threads = _num_threads;
}

template<class Engine, class Sseq> void MTPRNG<Engine,false,Sseq>::reseed(Sseq & _seq)
{
    std::mt19937 seeder(_seq);
    uint32_t per_engine_seed[seed_length];
    for(unsigned int j=0;j<seed_length;++j)
    {
        per_engine_seed[j] = seeder();
    }
    std::seed_seq derived_seed_seq(per_engine_seed, per_engine_seed+seed_length);
    engine.seed(derived_seed_seq);
};



//generic Sampler. All other sampler are derived from it.

template<class ET,bool MT, class Engine, class Sseq> //Sseq is supposed to satisfy the C++ concept "SeedSequence". The standard library has std::seed_seq as a canonical example.
                                                     //Engine is supposed to satisfy the C++ concept of a "Random Number Engine". <random> provides several of those, e.g. std::mt19937_64.
class Sampler
//Note :    In multi-threaded environment, we only have 1 sampler object. thread-number is given to sample();
{
    public:
    friend ostream & operator<< <ET,MT>(ostream &os, Sampler<ET,MT,Engine, Sseq>* const samplerptr);
    friend istream & operator>> <ET,MT>(istream &is, Sampler<ET,MT,Engine, Sseq>* const samplerptr);

    Sampler<ET,MT,Engine,Sseq> (Sseq & initial_seed): engine(initial_seed), sieveptr(nullptr)                      {}
    //We call init first, then custom_init (via init).
    void init(Sieve<ET,MT> * const sieve);
    virtual ~Sampler()=0; //needs to be virtual
    virtual SamplerType  sampler_type() const {return SamplerType::user_defined;};    //run-time type information.
                                                                    //This may be used to determine how to interpret a dump file.
                                                                    //defaults to user-defined.
                                                                    //Other values mean that the GaussSieve dumping routine is aware of the type, simplifying the syntax for dumping / reading.
    //TODO : Write some explanation how to do that.
    virtual LatticePoint<ET> sample(int thread=0)=0; //thread is the index of the calling thread (we need to keep separate PRNGs for each thread)
    //TODO : Allow sampling in subspaces, updating basis.

    private:
    virtual void custom_init()                                                                  {}         //called before any points are sampled;
    virtual ostream & dump_to_stream(ostream &os)  {return os;};    //dummy implementation of << operator.
    virtual istream & read_from_stream(istream &is){return is;};    //dummy implementation of >> operator.
    protected:
    MTPRNG<Engine, MT, Sseq> engine; //or engines
    Sieve<ET,MT> * sieveptr; //pointer to parent sieve.
};

template <class ET,bool MT, class Engine, class Sseq> Sampler<ET,MT, Engine,Sseq>::~Sampler() {} //actually needed, even though destructor is pure virtual as the base class destructor is eventually called implicitly.

template <class ET,bool MT, class Engine, class Sseq> void Sampler<ET,MT,Engine,Sseq>::init(Sieve<ET,MT> * const sieve)
{
    sieveptr = sieve;
    engine.init(sieve->get_num_threads());
    custom_init();
}

template<class ET,bool MT, class Engine, class Sseq> ostream & operator<<(ostream &os,Sampler<ET,MT,Engine,Sseq>* const samplerptr){return samplerptr->dump_to_stream(os);};
template<class ET,bool MT, class Engine, class Sseq> istream & operator>>(istream &is,Sampler<ET,MT,Engine,Sseq>* const samplerptr){return samplerptr->read_from_stream(is);};

/**
 * sampling integral Gaussians by rejection sampling
 */

 //Z must be an integral POD class

template<class Z, class Engine> Z GaussSieve::sample_z_gaussian(double s, double const center, Engine & engine, double const cutoff)
{
//Note : The following allows to access / modify floating point exceptions and modes.
//#pragma STDC FENV_ACCESS on
//This is too compiler/implementation-specific and does not work most of the time...

    static_assert(is_integral<Z>::value,"Return type for sample_z_gaussian must be POD integral type.");
    Z maxdev = static_cast<Z>(std::ceil(s * cutoff)); //maximum deviation of the Gaussian from the center. Note that maxdev may be 1.
    std::uniform_int_distribution<Z> uniform_in_range (std::floor(center-maxdev),std::ceil(center+maxdev));
    std::uniform_real_distribution<double> rejection_test(0.0,1.0); //defaults to value from [0,1), used in rejection sampling.
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

template<class Z, class Engine> Z GaussSieve::sample_z_gaussian_VMD(double const s2pi, double const center, Engine & engine, double const maxdeviation)
{
//Note : The following allows to access / modify floating point exceptions and modes.
//#pragma STDC FENV_ACCESS on
//This is too compiler/implementation-specific and does not work most of the time...

    static_assert(is_integral<Z>::value,"Return type for sample_z_gaussian must be POD integral type.");
    std::uniform_int_distribution<Z> uniform_in_range (std::floor(center-maxdeviation),std::ceil(center+maxdeviation));
    std::uniform_real_distribution<double> rejection_test(0.0,1.0); //defaults to value from [0,1), used in rejection sampling.
    Z closest_int = std::round(center); //closest int to center, i.e. most likely value.
    double adj = -(center-closest_int)*(center-closest_int); //negative squared distance to most likely value. Used to scale up the Gaussian weight function s.t. it is 1 at the most likely value.

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

        if( rejection_test(engine) <  std::exp(-std::fma(dist,dist,adj)/s2pi))
        {
            //std::feclearexcept(FE_UNDERFLOW);
            //std::feupdateenv(&env);
            return result;
        }
    }
}

#endif
