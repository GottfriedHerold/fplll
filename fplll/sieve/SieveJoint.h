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


#ifdef DO_INCLUDE_SIEVE_JOINT_H


template<class T>
class IsLongerVector_class;
template<class T>
class TerminationConditions;
template<class ET, bool MultiThreaded>
class Sieve;

#ifndef SIEVE_JOINT_H
//The following should be done only once (even if we include this header twice)
#include "LatticePoint.h"
#include "PointList.h"
#include <iostream>
#include <type_traits>
#include <sys/stat.h>
#include <fstream>

template<class ET>
using SieveST = Sieve<ET,false>;

template<class ET>
using SieveMT = Sieve<ET,true>;

template<class ET>
class TerminationConditions
{
    public:
    TerminationConditions() : default_condition(true){}; //if this is set, we are to ignore all other values.
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
    bool do_we_use_default_condition() const {return default_condition;};
    unsigned long int get_allowed_collisions() const {return allowed_collisions_;};
    void set_allowed_collisions(unsigned long const &colls) {allowed_collisions_=colls;do_we_check_collisions_=true;return;};
    unsigned long int get_allowed_list_size() const {return allowed_list_size_;};
    void set_allowed_list_size(unsigned long const &maxsize){allowed_list_size_=maxsize;do_we_check_list_size_=true;return;};
    ET get_target_length() const {return target_length_;};  //in case ET = mpz, it won't work
    void set_target_length(ET const &new_target_length) {target_length_=new_target_length;do_we_check_length_=true;return;};
    private:
    bool do_we_check_collisions_;
    bool do_we_check_length_;
    bool do_we_check_list_size_;
    bool default_condition;
    unsigned long int allowed_collisions_;
    unsigned long int allowed_list_size_;
    ET target_length_;
};
#endif
//The following may be included once or twice (with different values for GAUSS_SIEVE_IS_MULTI_THREADED)


template<class ET> //ET : underlying entries of the vectors. Should be a Z_NR<foo> - type. Consider making argument template itself.
class Sieve<ET, GAUSS_SIEVE_IS_MULTI_THREADED >
{

//data types: Note that these may change, so use these typedef's rather than constructions with ET directly.

using LPType           = LatticePoint<ET>;
using MainQueueType    = std::priority_queue< LPType, std::vector<LPType>, IsLongerVector_class<ET> >;
using MainListType     = PointListSingleThreaded<ET>;
//using MainListType     = std::list<LPType>;
using LatticeBasisType = ZZ_mat<typename ET::underlying_data_type>;
using SamplerType      = KleinSampler<typename ET::underlying_data_type, FP_NR<double>> *; //TODO : Should be a class with overloaded operator() or with a sample() - member.;

public:
//constructors:
Sieve() = delete;
Sieve(Sieve const &old ) = delete;
Sieve(Sieve &&old) = default;
Sieve & operator=(Sieve const & old)=delete;
Sieve & operator=(Sieve &&old) = default; //movable, but not copyable.

Sieve(LatticeBasisType B, unsigned int k=2, TerminationConditions<ET> termcond = {}, unsigned int verbosity_=1, int seed_sampler = 0);
//TODO: dump_status_to_stream
//TODO: read_status_from_stream -> Make constructor

//destructor:
~Sieve() {delete sampler;};


#ifdef GAUSS_SIEVE_SINGLE_THREADED
static bool const class_multithreaded = false;
#else
static bool const class_multithreaded = true;
#endif //class_multithreaded is for introspection, is_multithreaded is what the caller wants (may differ if we dump and re-read with different params)

void run_2_sieve(); //actually runs the Gauss Sieve.
LPType get_SVP(); //obtains Shortest vector and it's length. If sieve has not yet run, start it.
void run(); //runs the sieve specified by the parameters.
void print_status(int verb = -1, std::ostream &out = cout) const;
 //prints status to out. verb overrides the verbosity unless set to -1.
void dump_status_to_stream(std::string const &outfilename, bool overwrite = false);



//getter / setter functions

int get_verbosity() const {return verbosity;};
void set_verbosity(int new_verbosity) {verbosity=new_verbosity;return;};
unsigned int get_lattice_rank() const {return lattice_rank;};
unsigned int get_ambient_dimension() const {return ambient_dimension;};
unsigned int get_k() const {return sieve_k;};
void set_k(unsigned int new_k) {sieve_k=new_k;return;};
bool is_multithreaded_wanted() const {return multi_threaded_wanted;}; //Note: No setter
LPType get_shortest_vector_found() const {return shortest_vector_found;};
ET get_best_length2() const {return get_shortest_vector_found().norm2; } //in case ET = mpz, it won't work -- ET is supposed to be copy-constructible, so Z_NR<mpz_t> should work.
bool check_whether_sieve_is_running() const {return (sieve_status==SieveStatus::sieve_status_running);};
unsigned long int get_number_of_collisions() const {return number_of_collisions;};
unsigned long int get_number_of_points_sampled() const {return number_of_points_sampled;};
unsigned long int get_number_of_points_constructed() const {return number_of_points_constructed;};

private:

//Use termination Condition to check whether we are done, based on statistics so far.
bool check_if_done();

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

public: TerminationConditions<ET> term_cond; private: //to avoid complicated (due to template hack) friend - declaration.

//results

enum class SieveStatus
{
  sieve_status_error  =  -1, //indicates an error (add error codes as neccessary)
  sieve_status_init   =  1, //we have initialized data (and may yet initialize some more, but sieve has not started
  sieve_status_running=  2, //sieve is currently running
  sieve_status_finished=100 //sieve has finished
} sieve_status;
LPType shortest_vector_found; //including its length

//statistics

unsigned long int number_of_collisions;
unsigned long int number_of_points_sampled;
unsigned long int number_of_points_constructed; //sampling  + succesful pairs
//length of shortest vector contained in shortest_vector_found

//TODO: total time spent?


};

//MOVE TO CPP FROM HERE:


template<class ET> //ET : underlying entries of the vectors. Should be a Z_NR<foo> - type. Consider making argument template itself.
Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::Sieve(LatticeBasisType B, unsigned int k, TerminationConditions<ET> termcond , unsigned int verbosity_, int seed_sampler): //move to cpp //TODO:MT
main_list(),
main_queue(),
original_basis(B),
lattice_rank(B.get_rows()),
ambient_dimension(B.get_cols()), //Note : this means that rows of B form the basis.
multi_threaded_wanted(GAUSS_SIEVE_IS_MULTI_THREADED),
sieve_k(k),
sampler(nullptr),
verbosity(verbosity_),
term_cond(termcond),
sieve_status(SieveStatus::sieve_status_init),
shortest_vector_found(), //TODO : initialize to meaningful value, e.g. any vector of B.
number_of_collisions(0),
number_of_points_sampled(0),
number_of_points_constructed(0)
{
 sampler = new KleinSampler<typename ET::underlying_data_type , FP_NR<double>>(B, verbosity, seed_sampler);
 //TODO : initialize term_condition to some meaningful default.
};

template<class ET>
void Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::dump_status_to_stream(std::string const &outfilename, bool overwrite)
{
if(!overwrite)
{
  //checks if file exists
  struct stat buffer;
  if (::stat(outfilename.c_str(), &buffer) == 0)
  {
  cerr << "Trying to dump to existing file without overwrite flag set. Aborting dump." << endl;
  return;
  }
  std::ofstream of(outfilename,std::ofstream::out || std::ofstream::trunc);
}
}



#define SIEVE_JOINT_H
#endif
