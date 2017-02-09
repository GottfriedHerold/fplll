#ifndef SIEVE_SINGLE_THREADED_H
#define SIEVE_SINGLE_THREADED_H

#include "LatticePoint.h"
#include "PointList.h"

template<class T>
class IsLongerVector_class;
template<class T>
class TerminationConditions;
template<class ZT>
class SieveST;


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

template<class ZT> //T : underlying entries of the vectors. Should be a Z_NR<foo> - type.
class SieveST
{
using LPType           = LatticePoint<ZT>;
using MainQueueType    = std::priority_queue< LPType, std::vector<LPType>, IsLongerVector_class<ZT> >;
using MainListType     = std::list<LPType>;
using LatticeBasisType = std::list<LPType>;
using SamplerType      = void*; //TODO : Should be a class with overloaded operator() or with a sample() - member.;
public:

private:
MainListType MainList;
MainQueueType MainQueue;
LatticeBasisType OriginalBasis;
int LatticeRank;
int AmbientDimension; //consider merging theses into a latticespec class.

TerminationConditions<ZT> term_cond;
bool check_if_done();

SamplerType sampler;
};
 //single-threaded sieves

#endif
