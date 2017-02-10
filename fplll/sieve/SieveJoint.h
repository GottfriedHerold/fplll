#undef DO_INCLUDE_SIEVE_JOINT_H
#ifndef GAUSS_SIEVE_IS_MULTI_THREADED
#error wrong usage of SieveJoint.h -- 1
#endif

#if GAUSS_SIEVE_IS_MULTI_THREADED == false

#if !defined(GAUSS_SIEVE_SINGLE_THREADED) || defined(GAUSS_SIEVE_MULTI_THREADED)
#error wrong usage of SieveJoint.h -- 2
#endif

#ifndef SIEVE_JOINT_H_ST
#define SIEVE_JOINT_H_ST
#define DO_INCLUDE_SIEVE_JOINT_H
#endif

#elif GAUSS_SIEVE_IS_MULTI_THREADED == true
#if defined(GAUSS_SIEVE_SINGLE_THREADED) || !defined(GAUSS_SIEVE_MULTI_THREADED)
#error wrong usage of SieveJoint.h -- 3
#endif

#ifndef SIEVE_JOINT_H_MT
#define SIEVE_JOINT_H_MT
#define DO_INCLUDE_SIEVE_JOINT_H
#endif

#endif

#include "LatticePoint.h"
#include "PointList.h"

#ifdef DO_INCLUDE_SIEVE_JOINT_H


#include "LatticePoint.h"
#include "PointList.h"

template<class T>
class IsLongerVector_class;
template<class T>
class TerminationConditions;
template<class ZT, bool MultiThreaded>
class Sieve;

#ifndef SIEVE_JOINT_H
//do this only once.
template<class ZT>
using SieveST = Sieve<ZT,false>;

template<class ZT>
using SieveMT = Sieve<ZT,true>;

template<class ZT>
class TerminationConditions
{
    public:
    TerminationConditions() = default;
    TerminationConditions(TerminationConditions const &old)=default;
    TerminationConditions(TerminationConditions && old)=default;
    TerminationConditions & operator=(TerminationConditions const &old)=default;
    TerminationConditions & operator=(TerminationConditions && old)=default;
    ~TerminationConditions()= default;
    enum class TerminationType //different basic types of termination check
      {
      CheckCollisions=1,
      CheckLength=2,
      CheckListSize=3
      };
    bool do_we_check_collisions() const {return do_we_check_collisions_;};
    bool do_we_check_lenght() const {return do_we_check_length_;};
    bool do_we_check_list_size() const {return do_we_check_list_size_;};
    unsigned long int get_allowed_collisions() const {return allowed_collisions_;};
    void set_allowed_collisions(unsigned long const &colls) {allowed_collisions_=colls;do_we_check_collisions_=true;return;};
    unsigned long int get_allowed_list_size() const {return allowed_list_size_;};
    void set_allowed_list_size(unsigned long const &maxsize){allowed_list_size_=maxsize;do_we_check_list_size_=true;return;};
    ZT get_target_length() const {return target_length_;};
    void set_target_length(ZT const &new_target_length) {target_length_=new_target_length;do_we_check_length_=true;return;};
    private:
    bool do_we_check_collisions_;
    bool do_we_check_length_;
    bool do_we_check_list_size_;
    unsigned long int allowed_collisions_;
    unsigned long int allowed_list_size_;
    ZT target_length_;
};
#endif

template<class ZT> //ZT : underlying entries of the vectors. Should be a Z_NR<foo> - type. Consider making argument template itself.
class Sieve<ZT, GAUSS_SIEVE_IS_MULTI_THREADED >
{
using LPType           = LatticePoint<ZT>;
using MainQueueType    = std::priority_queue< LPType, std::vector<LPType>, IsLongerVector_class<ZT> >;
using MainListType     = std::list<LPType>;
using LatticeBasisType = std::list<LPType>;
using SamplerType      = void*; //TODO : Should be a class with overloaded operator() or with a sample() - member.;
public:

Sieve() = default;
Sieve(Sieve const &old ) = delete;
Sieve(Sieve &&old) = default;
Sieve & operator=(Sieve const & old)=delete;
Sieve & operator=(Sieve &&old) = default; //movable, but not copyable.
~Sieve()=default;
//Sieve( ) //TODO : Construct from LatticeBasis and Term. Conditions.
//TODO: dump_status_to_stream
//TODO: read_status_from_stream -> Make constructor

#ifdef GAUSS_SIEVE_SINGLE_THREADED
static bool const class_multithreaded = false;
#else
static bool const class_multithreaded = true;
#endif //class_multithreaded is for introspection, is_multithreaded is what the caller wants (may differ if we dump and re-read with different params)

void run_2_sieve(); //actually runs the Gauss Sieve.
LPType get_SVP(); //obtains Shortest vector and it's length. If sieve has not yet run, start it.
void run(); //runs the sieve specified by the parameters.
void print_status(int verb = -1) const; //prints status to cout. verb override the verbosity unless set to -1.

//getter / setter functions

int get_verbosity() const {return verbosity;};
void set_verbosity(int new_verbosity) {verbosity=new_verbosity;return;};
unsigned int get_lattice_rank() const {return lattice_rank;};
unsigned int get_ambient_dimension() const {return ambient_dimension;};
unsigned int get_k() const {return sieve_k;};
void set_k(unsigned int new_k) {sieve_k=new_k;return;};
bool is_multithreaded_wanted() const {return multi_threaded_wanted;}; //Note: No setter
LPType get_shortest_vector_found() const {return shortest_vector_found;};
ZT get_best_length2() const {return get_shortest_vector_found().norm2; }
bool check_whether_sieve_is_running() const {return sieve_is_running;};
unsigned long int get_number_of_collisions() const {return number_of_collisions;};
unsigned long int get_number_of_points_sampled() const {return number_of_points_sampled;};
unsigned long int get_number_of_points_constructed() const {return number_of_points_constructed;};


private:

//Note: The member fields of Sieve denote the (global) "internal" status of the sieve during a run or execution.
//It should be possible to dump the status to harddisk and resume from dump using that information.
//It should also be possible to suspend the run of the sieve, change (certain) parameters (like k!) and resume.

//main data that is changing.

MainListType main_list;
MainQueueType main_queue;

//information about lattice and algorithm we are using

LatticeBasisType original_basis;
unsigned int lattice_rank;
unsigned int ambient_dimension; //consider merging theses into a latticespec class.
bool multi_threaded_wanted;
//unsigned int num_threads_wanted;
unsigned int sieve_k; //parameter k of the sieve currently running.
SamplerType sampler;
int verbosity;
public:
TerminationConditions<ZT> term_cond;
private:
//results

bool check_if_done(); //Use termination Condition to check whether we are done, based on statistics so far.
bool sieve_is_running;
LPType shortest_vector_found; //including its length

//statistics

unsigned long int number_of_collisions;
unsigned long int number_of_points_sampled;
unsigned long int number_of_points_constructed; //sampling  + succesful pairs
//length of shortest vector contained in shortest_vector_found

//TODO: total time spent?


};

#define SIEVE_JOINT_H
#endif
