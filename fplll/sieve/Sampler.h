#ifndef SAMPLER_H
#define SAMPLER_H

//forward declarations

template<class ET, bool MT>
class Sampler;
template<class ET,bool MT> ostream & operator<<(ostream &os, Sampler<ET,MT>* const samplerptr); //printing
template<class ET,bool MT> istream & operator>>(istream &is, Sampler<ET,MT>* const samplerptr); //reading (may also be used by constructor from istream)
enum class SamplerType
{
    user_defined = 0,
};
//includes


template<class ET,bool MT>
class Sampler
//Note :    In multi-threaded environment, we only have 1 sampler object
//          caller_thread is set to -1 if we make a single-threaded call
{
    public:
    friend ostream & operator<< <ET,MT>(ostream &os, Sampler<ET,MT>* const samplerptr);
    friend istream & operator>> <ET,MT>(istream &is, Sampler<ET,MT>* const samplerptr);
    virtual void init(Sieve<ET,MT> * const sieve, long int const seed) {}                           //called before any points are sampled;
    virtual ~Sampler()=0; //needs to be virtual
    virtual SamplerType  sampler_type() const {return SamplerType::user_defined;};    //run-time type information.
                                                                    //This is used to determine how to interpret a dump file.
                                                                    //defaults to user-defined.
                                                                    //Other values mean that the GaussSieve dumping routine is aware of the type, simplifying the syntax for dumping / reading.
    //TODO : Write some explanation how to do that.
    virtual LatticePoint<ET> sample(int thread=0); //thread is the index of the calling thread (we need to keep separate PRNGs for each thread)
    //TODO : Allow sampling in subspaces, updating basis.

    private:
    virtual ostream & dump_to_stream(ostream &os)  {return os;};    //dummy implementation of << operator.
    virtual istream & read_from_stream(istream &is){return is;};    //dummy implementation of >> operator.

};
template <class ET,bool MT>
Sampler<ET,MT>::~Sampler() {} //actually needed, even though destructor is pure virtual as the base class destructor is eventually called implicitly.

template<class ET,bool MT> ostream & operator<<(ostream &os,Sampler<ET,MT>* const samplerptr){return samplerptr->dump_to_stream(os);};
template<class ET,bool MT> istream & operator>>(istream &is,Sampler<ET,MT>* const samplerptr){return samplerptr->read_from_stream(is);};




#endif
