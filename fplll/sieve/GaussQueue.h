#ifndef GAUSS_QUEUE_H
#define GAUSS_QUEUE_H

/* defines the classes used for the main Queues in the Gauss Sieve */

template <class ET, bool MT>
class GaussQueue;

template <class ET>
using GaussQueueST=GaussQueue<ET,false>;

template <class ET>
using GaussQueueMT=GaussQueue<ET,true>;

#include "LatticePoint.h"
#include <mutex>
#include <atomic>
#include <queue>
#include "assert.h"
#include "SieveGauss.h"

template<class ET> class IsLongerVector_classPtr
{
    public: bool operator() (LatticePoint<ET>* const &A, LatticePoint<ET>* const & B)
    {
     return (A->get_norm2() > B->get_norm2() );
    }
};

template<class ET> //single-threaded version:
class GaussQueue<ET,false>
{
public:
    using LPType = LatticePoint<ET>;
    using QueueType =      std::priority_queue< LPType* , std::vector<LPType* >, IsLongerVector_classPtr<ET> >;
    using SamplerType =    KleinSampler<typename ET::underlying_data_type, FP_NR<double> > ;


    GaussQueue()=delete;
    GaussQueue(Sieve<ET,false> *caller_sieve);
    GaussQueue(GaussQueue const &old) = delete;
    GaussQueue(GaussQueue &&old) = delete;
    GaussQueue& operator= (GaussQueue const &old)=delete;
    GaussQueue& operator= (GaussQueue &&old) = delete;
    ~GaussQueue();

    void pop();
    LPType const & top() const;
    bool empty() const; //we might as well always return false!
    typename QueueType::size_type size() const;


private:
    QueueType our_queue;
    Sieve<ET,false>* gauss_sieve; //caller object.
    SamplerType *sampler; //controlled by the GaussSieve currently?
};

template<class ET>
void GaussQueue<ET,false>::pop()
{
    if(!our_queue.empty())
    {
        delete our_queue.front();
        our_queue.pop();
    }
}

template<class ET>
typename GaussQueue<ET,false>::LPType const & GaussQueue<ET,false>::top() const
{
    if(our_queue.empty())
    {
        LPType * newpointptr = new LPType ( sampler->sample() );
    }
    return * ( our_queue.top() );
}

template<class ET>
typename GaussQueue<ET,false>::QueueType::size_type GaussQueue<ET,false>::size() const
{
    return our_queue.size();
}

template<class ET>
bool GaussQueue<ET,false>::empty() const
{
    return our_queue.empty();
}

template<class ET>
GaussQueue<ET,false>::~GaussQueue()
{
//TODO: Delete sampler if owned.
//Note that deletion of points is already done.
//May need to delete points if underlying type is switched to pointer.
}

#endif
