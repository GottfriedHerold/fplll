
/*
------DO NOT INCLUDE THIS FILE MANUALLY.--------
USE SieveGauss.h / SieveMT.h / SieveST.h INSTEAD

*/

/* Preprocessor magic starts here */
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

//
//end of (most) preprocessor stuff
//




#ifndef SIEVE_JOINT_H
/*

THE FOLLOWING PARTS ARE ONLY PARSED ONCE BY THE COMPILER,
EVEN IF WE USE BOTH MULTI-THREADED AND SINGLE-THREADED VARIANTS.

DECLARATIONS CAN GO HERE.

HELPER CLASSES AND FUNCTIONS WHICH ARE NOT TEMPLATES WITH TEMPLATE PARAMETER
GAUSS_SIEVE_IS_MULTI_THREADED
NEED TO GO HERE:
*/

/*INCLUDES */

#include "LatticePoint.h"
#include "PointList.h"
#include <iostream>
#include <type_traits>
#include <sys/stat.h>
#include <fstream>
#include <string>
#include <exception>

/*Declarations */

template<class T>
class IsLongerVector_class;
template<class T>
class TerminationConditions;
template<class ET, bool MultiThreaded>
class Sieve;


bool string_consume(istream &is, std::string const & str, bool elim_ws= true, bool verbose=true);

template<class ET>
using SieveST = Sieve<ET,false>;

template<class ET>
using SieveMT = Sieve<ET,true>;


/* CLASS FOR TERMINATION CONDITIONS */

class bad_dumpread_TermCond:public std::runtime_error{public:bad_dumpread_TermCond():runtime_error("Dump read failed for Termination Condition"){}}; //exception indicating that read from dump failed.
template<class ET> ostream & operator<<(ostream &os,TerminationConditions<ET> const &term_cond); //printing
template<class ET> istream & operator>>(istream &is,TerminationConditions<ET> &term_cond); //reading (also used by constructor from istream)
template<class ET>
class TerminationConditions
{
    friend ostream & operator<< <ET>(ostream &os,TerminationConditions const &term_cond); //printing
    friend istream & operator>> <ET>(istream &is,TerminationConditions &term_cond);
    public:
    TerminationConditions() : default_condition(true){}; //if this is set, we are to ignore all other values.
    explicit TerminationConditions(istream &is):default_condition(true){is>> (*this);};
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
}; //end of termination condition class


#endif

/*
EVERYTHING BELOW HERE IS POTENTIALLY INCLUDED TWICE.
TEMPLATES WITH TEMPLATE ARGUMENT GAUSS_SIEVE_IS_MULTI_THREADED
GO HERE.
*/

//The following may be included once or twice (with different values for GAUSS_SIEVE_IS_MULTI_THREADED)

template<class ET>
class Sieve<ET, GAUSS_SIEVE_IS_MULTI_THREADED >
{
/*DATA TYPES*/
using LPType           = LatticePoint<ET>;
using MainQueueType    = std::priority_queue< LPType, std::vector<LPType>, IsLongerVector_class<ET> >;
using MainListType     = PointListSingleThreaded<ET>;
//using MainListType     = std::list<LPType>;
using LatticeBasisType = ZZ_mat<typename ET::underlying_data_type>;
using SamplerType      = KleinSampler<typename ET::underlying_data_type, FP_NR<double>> *; //TODO : Should be a class with overloaded operator() or with a sample() - member.;

public:
/*CONSTRUCTORS / DESTRUCTORS */
Sieve() = delete;
Sieve(Sieve const &old ) = delete;
Sieve(Sieve &&old) = default;
Sieve & operator=(Sieve const & old)=delete;
Sieve & operator=(Sieve &&old) = default; //movable, but not copyable.
explicit Sieve(LatticeBasisType B, unsigned int k=2, TerminationConditions<ET> termcond = {}, unsigned int verbosity_=1, int seed_sampler = 0);
//explicit Sieve(std::string const &infilename); //read from dump file.
~Sieve() {delete sampler;};

#ifdef GAUSS_SIEVE_SINGLE_THREADED
static bool const class_multithreaded = false;
#else
static bool const class_multithreaded = true;
#endif //class_multithreaded is for introspection, is_multithreaded is what the caller wants (may differ if we dump and re-read with different params)

    void run_2_sieve(); //actually runs the Gauss Sieve.
    void SieveIteration (LatticePoint<ET> &p); //one run through the main_list
    
LPType get_SVP(); //obtains Shortest vector and it's length. If sieve has not yet run, start it.
void run(); //runs the sieve specified by the parameters.
//void print_status(int verb = -1, std::ostream &out = cout) const{dump_status_to_stream(out,false,verb);};
 //prints status to out. verb overrides the verbosity unless set to -1.
void dump_status_to_file(std::string const &outfilename, bool overwrite = false); //dumps to file
void dump_status_to_stream(ostream &of, bool everything = false, int verb=-1); //dumps to stream. Can be read back if everything == true. Otherwise, verbosity determines what is output.


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
unsigned long int get_current_list_size() const{return current_list_size;};
unsigned long int get_current_queue_size()const{return main_queue.size();};

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
unsigned int ambient_dimension; //consider merging these into a latticespec struct.
bool multi_threaded_wanted;
//unsigned int num_threads_wanted;
unsigned int sieve_k; //parameter k of the sieve currently running.
SamplerType sampler;
int verbosity;

public: TerminationConditions<ET> term_cond; private: //to avoid complicated (due to template hack) friend - declaration.

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
unsigned long int current_list_size;
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
number_of_points_constructed(0),
current_list_size(0)
{
 sampler = new KleinSampler<typename ET::underlying_data_type , FP_NR<double>>(B, verbosity, seed_sampler);
 //TODO : initialize term_condition to some meaningful default.
};



/*DUMPING / READING ROUTINES */

//Note: Actually, we want an unformatted binary dump. Unfortunately, the underlying FPLLL types support only
// << and >> operations with formated input / output. So we will do with formatted input / output for now.
//The main issue here is that mpz_t only provides formatted stream - I/O or unformatted I/O via old-style C FILE* interfaces, neither of which is what we really want.

// We assume that << writing data to a filestream and >> reading it back will give back the same data.

// This might cause problems if e.g.:
// separation chars (whitespace) are used in a data field
// locales are different
// some subobject changes format flags
// output loses data (e.g. rounding of floats)

//helper functions:

/*
Reads length(str) chars from stream is, expecting them to equal str. If what is read differs we output false. If verbose, we also display an error.
*/

#ifndef SIEVE_JOINT_H

bool string_consume(istream &is, std::string const & str, bool elim_ws, bool verbose)
{
  unsigned int len = str.length();
  char *buf = new char[len+1];
  buf[len] = 0; //for error message.
  if (elim_ws)
  {
    is >> std::ws;
  }
  is.read(buf,len);
  if(is.gcount() != len)
  {
    if(verbose)
    {
      cerr << "Failure reading header: Expected to read" << str << endl;
      cerr << "Read only "<<is.gcount() << "bytes. String read was" << buf<<endl;
    }
    return false;
  }
  if(elim_ws)
  {
    is >> std::ws;
  }
  if(str.compare(0,len,buf,len)!=0)
  {
    if(verbose)
    {
      cerr << "Failure reading header: Expected to read" << str << endl;
      cerr << "Read instead:" << buf << endl;
    }
    return false;
  }
  return true;
}

template<class ET> ostream & operator<<(ostream &os,TerminationConditions<ET> const &term_cond) //printing
{

  os << "Default_Conditions=" << term_cond.default_condition << endl;
  if(!term_cond.default_condition)
  {
    os << "Check Collisions=" << term_cond.do_we_check_collisions_ << endl;
    if(term_cond.do_we_check_collisions_)
    {
      os << "Number=" << term_cond.allowed_collisions_ << endl;
    }
    os << "Check List Size=" << term_cond.do_we_check_list_size_ << endl;
    if(term_cond.do_we_check_list_size_)
    {
      os << "Number=" << term_cond.allowed_list_size_ << endl;
    }
    os << "Check Target Length=" << term_cond.do_we_check_length_ << endl;
    if(term_cond.do_we_check_length_)
    {
      os << "Target Length=" << term_cond.target_length_ << endl;
    }
  }
return os;
}
template<class ET> istream & operator>>(istream &is,TerminationConditions<ET> &term_cond)
{
    bool do_we_check_collisions_;
    bool do_we_check_length_;
    bool do_we_check_list_size_;
    bool default_condition;
    unsigned long int allowed_collisions_;
    unsigned long int allowed_list_size_;
    ET target_length_;
 //We should probably throw an exception rather than return is.
  if (!string_consume(is,"Default_Conditions=")) throw bad_dumpread_TermCond();
  is >> term_cond.default_condition;
  if(!term_cond.default_condition)
  {
    if(!string_consume(is,"Check Collisions=")) throw bad_dumpread_TermCond();
    is>> term_cond.do_we_check_collisions_;
    if(term_cond.do_we_check_collisions_)
    {
      if(!string_consume(is,"Number=")) throw bad_dumpread_TermCond();
      is >> term_cond.allowed_collisions;
    }
    if(!string_consume(is,"Check List Size=")) throw bad_dumpread_TermCond();
    is>> term_cond.do_we_check_list_size_;
    if(term_cond.do_we_check_list_size_)
    {
      if(!string_consume(is,"Number=")) throw bad_dumpread_TermCond();
      is>> term_cond.allows_list_size_;
    }
    if(!string_consume(is,"Check Target Length=")) throw bad_dumpread_TermCond();
    is>>term_cond.do_we_check_length_;
    if(term_cond.do_we_check_length_)
    {

    }

  }
return is;

} //reading (also used by constructor from istream)


#endif


#define SIEVE_FILE_ID "kTuple-Sieve dump file"
//version string for dump file
#define SIEVE_VER_STR "Version TEST1"
#define SIEVE_FILE_ML "Main List"
#define SIEVE_FILE_QUEUE "Queue"


template<class ET>
void Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::dump_status_to_file(std::string const &outfilename, bool overwrite)
{
if(verbosity>=2)
{
  cout << "Dumping to file " << outfilename << " ..." << endl;
}
if(!overwrite)
{
  //checks if file exists
  struct stat buffer;
  if (::stat(outfilename.c_str(), &buffer) == 0)
  {
  if (verbosity>=1)
  {
    cerr << "Trying to dump to existing file without overwrite flag set. Aborting dump." << endl;
  }
  return;
  }
}
std::ofstream of(outfilename,std::ofstream::out | std::ofstream::trunc);
dump_status_to_stream(of, true);
if(verbosity>=2)
{
  cout << "Dump successful." << endl;
}
}

template<class ET>
void Sieve<ET,GAUSS_SIEVE_IS_MULTI_THREADED>::dump_status_to_stream(ostream &of, bool everything, int verb)
{
int howverb = everything ? 100 : (verb==-1 ? verbosity : verb);
if(howverb>=2) of << SIEVE_FILE_ID << endl;
if(howverb>=2) of << SIEVE_VER_STR << endl;
if(howverb>=2) of << "--Params--" << endl;
if(howverb>=2) of << "Multithreaded version=" << class_multithreaded << endl;
if(howverb>=2) of << "Multithreaded is wanted" << multi_threaded_wanted << endl;
if(howverb>=2) of << "k=" << sieve_k << endl;
if(howverb>=2) of << "verbosity=" << verbosity << endl;
if(howverb>=1) of << "sieve_status=" << static_cast<int>(sieve_status) << endl;
if(howverb>=2) of << "Lattice rank=" << lattice_rank << endl;
if(howverb>=2) of << "Ambient dimension=" << ambient_dimension << endl;
if(howverb>=2) of << "Termination Conditions:" << endl;
if(howverb>=2) of << term_cond;
if(howverb>=2) of << "Sampler:"<< endl;
if(howverb>=2) of << "Sampler Initialized" << static_cast<bool>(sampler!=nullptr) << endl;
if(sampler!=nullptr)
{
if(howverb>=2) cerr << "Note : Dumping of internal data of sampler not yet supported" << endl;
  //dump internals of sampler?
}
if(howverb>=1) of << "Original Basis:" << endl;
if(howverb>=1) of << original_basis;
if(howverb>=2) of << "--End of Params--" << endl << endl;
if(howverb>=1) of << "--Statistics--" << endl;
if(howverb>=1) of << "Number of collisions=" << number_of_collisions << endl;
if(howverb>=1) of << "Number of Points Sampled=" << number_of_points_sampled << endl;
if(howverb>=1) of << "Number of Points Constructed=" << number_of_points_constructed << endl;
if(howverb>=1) of << "Best vector found so far=" << shortest_vector_found << endl; //TODO : Display length seperately
if(howverb>=1) of << "Current List Size=" << get_current_list_size() << endl;
if(howverb>=1) of << "Current Queue Size="<< get_current_queue_size()<< endl;
if(howverb>=1) of << "--End of Statistics--" << endl << endl;
if(howverb>=3)
{
  of << "--Main List--" << endl;
  for(auto it = main_list.begin();it!=main_list.end();++it)
  {
    of << (*it);
  }
  of << "--End of Main List--" << endl << endl;
  of << "--Main Queue--" << endl;
//  for(auto it = main_queue.begin();it!=main_queue.end();++it)
  cerr << "Dumping of main queue not supported yet.";
  of << "--End of Main Queue--";
  {

  }
}


}

#define SIEVE_JOINT_H
#endif
