#ifndef TERM_COND_NEW_H
#define TERM_COND_NEW_H

//TODO : Change to this version later.
//Reconsider :
//pure virtual function version breaks dumping / reading (without some RTT)
//templated version locks into static typing. However, type of term-cond might depend on user input.

//solution : provide limited form of RTT.




//Note : Removed below idea, because we want run-time typing.
/* Concept "TerminationCondition":

A termination condition is a template<class ET> class TerminationCondition with at least the following public members: (see DummyTerminationCondition for an example)

- public member typedef simple, set to either std::true_type or std::false_type (Note : May need to #include <type_traits> )
- public default constructor
- public member function templates
template<bool MT> void init(GaussSieve<ET,MT> const & sieve);
template<bool MT> int check(GaussSieve<ET,MT> const & sieve);
template<bool MT> int check_vec(GaussSieve<ET,MT> const & sieve, ET const & norm2);

- non-member function templates
template<class ET> ostream & operator<<(ostream &os,TerminationCondition<ET> const &term_cond);
template<class ET> istream & operator>>(istream &is,TerminationCondition<ET> &term_cond);

It needs to have the following semantics:

init(...) is called if the GaussSieve is started. You may assume that the parameters do no longer change. Note that init() may be called multiple times if the sieve is suspended and parameters change.
check(...) returns 1 if sieve is considered finished, 0 otherwise. Needs to be reentrant if MT==true.
check_vec(...) has the same semantics as check, but is called whenever a new lattice point of norm^2 norm2 is found. Needs to be reentrant if MT==true.

If simple==true_type, we assume that check() is a function of length only. In this case, we only call check_vec() if a new shortest (so far) vector is found.
check() may be called rarely.

Note : Output of check functions is int rather than bool to enable future extensions to output return values indicating "suspend" or "dump"

The stream operators are for dumping / reading.
We assume that if we dump a Termination condition T1 to a filestream and read it into T2, then T2 will be in the same status as T1.
TODO: reading is unused so far.
*/
template<class ET, bool MT>
class TerminationCondition;
template<class ET,bool MT> ostream & operator<<(ostream &os,TerminationCondition<ET,MT>* const term_cond); //printing
template<class ET,bool MT> istream & operator>>(istream &is,TerminationCondition<ET,MT>* const term_cond); //reading (also used by constructor from istream)


#include "SieveGauss.h"

/*
    example dummy Termination condition below
*/

//template<class ET> ostream & operator<<(ostream &os,DummyTerminationCondition<ET> const &term_cond); //printing
//template<class ET> istream & operator>>(istream &is,DummyTerminationCondition<ET> &term_cond); //reading (also used by constructor from istream)
//template<class ET>
//class DummyTerminationCondition //This class serves only to exemplify the interface, lacking support for "Concepts" in any compiler execpt GCC (experimental) atm.
//{
//    public:
//    template<class ET> ostream & operator<<(ostream &os,TerminationConditions<ET> const &term_cond); //printing
//    template<class ET> istream & operator>>(istream &is,TerminationConditions<ET> &term_cond); //reading (also used by constructor from istream)
//    using simple = std::true_type;
//
//    DummyTerminationCondition() = default;
//
//    template<bool MT>
//    void init(GaussSieve<ET,MT> const & sieve) = delete;
//
//    template<bool MT>
//    int check(GaussSieve<ET,MT> const & sieve) = delete;
//
//    template<bool MT>
//    int check_vec(GaussSieve<ET,MT> const & sieve, ET const & length) = delete;
//}

template<class ET,bool MT>
class TerminationCondition //This class serves only to exemplify the interface, lacking support for "Concepts" in any compiler execpt GCC (experimental) atm.
{
    public:
    friend ostream & operator<< <ET,MT>(ostream &os,TerminationCondition<ET,MT>* const &term_cond);
    friend istream & operator>> <ET,MT>(istream &is,TerminationCondition<ET,MT>* &term_cond);
    TerminationCondition(): buf(nullptr){};
    virtual void init(Sieve<ET,MT> const & sieve) = 0;
    virtual int check(Sieve<ET,MT> const & sieve) = 0;
    virtual int check_vec(Sieve<ET,MT> const & sieve, ET const & length) = 0;

    virtual bool is_simple() const {return false;};
    virtual int  termination_condition_type() const {return 0;}; //run-time type information. This is used to determine how to interpret a dump file. 0 means user-defined (which is the default).
    //TODO : Write some explanation how to do that.
    virtual ~TerminationCondition() = 0;

    private:
    virtual ostream & dump_to_stream(ostream &os);    //implementation of << operator.
    virtual istream & read_from_stream(istream &is);  //implementation of >> operator.
    unsigned long data_size;        //holds user-defined data, needed to forward dumping / reading routines to caller.
    char * buf;

};
template <class ET,bool MT>
TerminationCondition<ET,MT>::~TerminationCondition() {if(buf!=nullptr)delete[] buf;} //actually needed, as the base class destructor is eventually called implicitly.

template<class ET,bool MT> ostream & operator<<(ostream &os,TerminationCondition<ET,MT>* const &term_cond){return term_cond->dump_to_stream(os);};
template<class ET,bool MT> istream & operator>>(istream &is,TerminationCondition<ET,MT>* &term_cond){return term_cond->read_from_stream(is);};


#endif
