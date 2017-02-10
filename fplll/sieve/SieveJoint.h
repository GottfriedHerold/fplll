#ifndef SIEVE_SINGLE_THREADED_H
#define SIEVE_SINGLE_THREADED_H

#include "LatticePoint.h"
#include "PointList.h"

template<class T>
class IsLongerVector_class;
template<class T>
class TerminationConditions;
template<class ZT, bool MultiThreaded>
class Sieve;

template<class ZT>
using SieveST = Sieve<ZT,false>;

template<class ZT>
using SieveMT = Sieve<ZT,true>;

template<class ZT>
class TerminationConditions
{
    int NumberOfCollisions;
    bool CheckCollisions;
    ZT TargetLength;
    bool CheckLength;
    int ListSizeLimit;
    bool CheckListSize;
};

template<class ZT> //ZT : underlying entries of the vectors. Should be a Z_NR<foo> - type. Consider making argument template itself.
class Sieve<ZT, false>
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

//TODO: dump_status_to_stream
//TODO: read_status_from_stream -> Make constructor

void run_2_sieve(); //actually runs the Gauss Sieve.
ZPType get_SVP(); //obtains Shortest vector and it's length. If sieve has not yet run, start it.
void run(); //runs the sieve specified by the parameters.

private:

//Note: The member fields of Sieve denote the (global) "internal" status of the sieve during a run or execution.
//It should be possible to dump the status to harddisk and resume from dump using that information.
//It should also be possible to suspend the run of the sieve, change (certain) parameters (like k!) and resume.

MainListType MainList;
MainQueueType MainQueue;
LatticeBasisType OriginalBasis;
int lattice_rank;
int ambient_dimension; //consider merging theses into a latticespec class.
bool multi_threaded;
int sieve_k; //parameter k of the sieve currently running.

TerminationConditions<ZT> term_cond;
bool check_if_done();
bool sieve_is_running;
LPType shortest_vector_found;

SamplerType sampler;
};

#endif
