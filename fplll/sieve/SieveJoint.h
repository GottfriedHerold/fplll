/*
------DO NOT INCLUDE THIS FILE MANUALLY.--------
------     USE SieveGauss.h INSTEAD     --------
*/

//SIEVE_JOINT_H_ST and SIEVE_JOINT_H_MT are separate include guards.
//They are set if we already included this file (for the ST / MT case) respectively.
//There is also an include guard SIEVE_JOINT_H, which is set _after_ we included this file at least once already.
//Use this to condition on the second pass.

#undef DO_INCLUDE_SIEVE_JOINT_H
#if GAUSS_SIEVE_IS_MULTI_THREADED == false

#ifndef SIEVE_JOINT_H_ST
#define SIEVE_JOINT_H_ST
#define DO_INCLUDE_SIEVE_JOINT_H
#endif

#elif GAUSS_SIEVE_IS_MULTI_THREADED == true
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

DECLARATIONS MAY GO HERE OR TO SieveGauss.h.

HELPER CLASSES AND FUNCTIONS WHICH ARE NOT TEMPLATES WITH TEMPLATE PARAMETER
GAUSS_SIEVE_IS_MULTI_THREADED
NEED TO GO HERE OR TO SieveGauss.h:
*/

/*
  Forward Declarations
  These go before includes to allow cyclic dependencies
*/

template<class T>
class IsLongerVector_class; //class wrapper to compare vectors by length

template<class ET, bool MultiThreaded>
class Sieve;

bool string_consume(istream &is, std::string const & str, bool elim_ws= true, bool verbose=true); //helper function for dumping/reading


/*INCLUDES */

#include "LatticePoint.h"
#include "PointList.h"
#include <iostream>
#include <type_traits>
#include <sys/stat.h>
#include <fstream>
#include <string>
#include <exception>
#include "TermCond.h"
#include "GaussQueue.h"

#endif //end of ONLY-ONCE part

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
    //using MainQueueType =std::queue<LPType>;
    //using MainQueueType = GaussQueue<ET,GAUSS_SIEVE_IS_MULTI_THREADED>;
    using MainListType     = PointListSingleThreaded<ET>;
    //using MainListType     = std::list<LPType>;
    using LatticeBasisType = ZZ_mat<typename ET::underlying_data_type>;
    using SamplerType      = KleinSampler<typename ET::underlying_data_type, FP_NR<double>> *; //TODO : Should be a class with overloaded operator() or with a sample() - member.;

public:
    /*FRIENDS */
    friend GaussQueue<ET,GAUSS_SIEVE_IS_MULTI_THREADED>;
    /*CONSTRUCTORS / DESTRUCTORS */
    Sieve() = delete;
    Sieve(Sieve const &old ) = delete;
    Sieve(Sieve &&old) = default;
    Sieve & operator=(Sieve const & old)=delete;
    Sieve & operator=(Sieve &&old) = default; //movable, but not copyable.
    explicit Sieve(LatticeBasisType B, unsigned int k=2, TerminationConditions<ET> termcond = {}, unsigned int verbosity_=1, int seed_sampler = 0);
    //explicit Sieve(std::string const &infilename); //read from dump file.
    ~Sieve()
    {
        delete sampler;
    };

#ifdef GAUSS_SIEVE_SINGLE_THREADED
    static bool const class_multithreaded = false;
#else
    static bool const class_multithreaded = true;
#endif //class_multithreaded is for introspection, is_multithreaded is what the caller wants (may differ if we dump and re-read with different params)

    void run_2_sieve(ET target_norm); //actually runs the Gauss Sieve.
    void SieveIteration (LatticePoint<ET> &p); //one run through the main_list
    LPType get_SVP(); //obtains Shortest vector and it's length. If sieve has not yet run, start it.
    void run(); //runs the sieve specified by the parameters.
    void print_status(int verb = -1, std::ostream &out = cout) {dump_status_to_stream(out,false,verb);};
    //prints status to out. verb overrides the verbosity unless set to -1.
    void dump_status_to_file(std::string const &outfilename, bool overwrite = false); //dumps to file
    void dump_status_to_stream(ostream &of, bool everything = false, int verb=-1); //dumps to stream. Can be read back if everything == true. Otherwise, verbosity determines what is output.

//getter / setter functions

    int get_verbosity() const                                   {return verbosity;};
    void set_verbosity(int new_verbosity)                       {verbosity=new_verbosity;return;};
    unsigned int get_lattice_rank() const                       {return lattice_rank;};
    unsigned int get_ambient_dimension() const                  {return ambient_dimension;};
    unsigned int get_k() const                                  {return sieve_k;};
    void set_k(unsigned int new_k)                              {sieve_k=new_k;return;};
    bool is_multithreaded_wanted() const                        {return multi_threaded_wanted;}; //Note: No setter
    LPType get_shortest_vector_found() const                    {return shortest_vector_found;};
    ET get_best_length2() const                                 {return get_shortest_vector_found().norm2;};
    bool check_whether_sieve_is_running() const                 {return (sieve_status==SieveStatus::sieve_status_running);};
    unsigned long int get_number_of_collisions() const          {return number_of_collisions;};
    unsigned long int get_number_of_points_sampled() const      {return number_of_points_sampled;};
    unsigned long int get_number_of_points_constructed() const  {return number_of_points_constructed;};
    unsigned long int get_current_list_size() const             {return current_list_size;};
    unsigned long int get_current_queue_size()const             {return main_queue.size();};

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

public:
    TerminationConditions<ET> term_cond;
private: //to avoid complicated (due to template hack) friend - declaration.

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

template<class ET>
bool check2red (LatticePoint<ET> &p1, const LatticePoint<ET> &p2);


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




#endif // SIEVE_JOINT_H


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
    if(howverb>=2) of << "Original Basis:" << endl;
    if(howverb>=2) of << original_basis;
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
        for(auto it = main_list.begin(); it!=main_list.end(); ++it)
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
#endif // DO_INCLUDE_SIEVE_JOINT_H
