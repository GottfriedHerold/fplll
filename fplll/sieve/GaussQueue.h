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
#include <utility>

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
    #ifndef USE_REGULAR_QUEUE
    using QueueType =      std::priority_queue< LPType* , std::vector<LPType* >, IsLongerVector_classPtr<ET> >;
    #else
    using QueueType =      std::queue<LPType*>
    #endif
    using SamplerType =    KleinSampler<typename ET::underlying_data_type, FP_NR<double> > ;
    GaussQueue()=delete;
    GaussQueue(Sieve<ET,false> *caller_sieve); //only constructor
    GaussQueue(GaussQueue const &old) = delete;
    GaussQueue(GaussQueue &&old) = delete;
    GaussQueue& operator= (GaussQueue const &old)=delete;
    GaussQueue& operator= (GaussQueue &&old) = delete;
    ~GaussQueue();

    bool empty() const; //we might as well always return false!
    typename QueueType::size_type size() const; //returns size of queue (used for diagnostics and statistics only)
    void push(LPType const &val); //puts a copy of val in the queue
    void push(LPType && val);     //uses move semantics for that.
    void give_ownership(LPType * const valptr); //takes a pointer to a list point and puts the point into the queue, moves ownership (avoids copying)
    LPType  true_pop(); //removes front element from queue *and returns it*.
    LPType* pop_take_ownership(); //removes front elements from queue and returns handle to it. Transfers ownership to the caller. Return type might change, but should be dereferencable, deleteable.

private:
    QueueType main_queue;
    Sieve<ET,false>* gauss_sieve; //caller object.
    SamplerType *sampler; //controlled by the GaussSieve currently. TODO: Change that
};

template<class ET>
GaussQueue<ET,false>::GaussQueue( Sieve<ET,false> *caller_sieve)  //constructor
:
main_queue(),
gauss_sieve(caller_sieve),
sampler(nullptr)
{
    assert(caller_sieve!=nullptr);
    sampler=caller_sieve->sampler; //TODO: Remove sampler from SieveJoint.h and place it under control of the queue.
}

template<class ET>
typename GaussQueue<ET,false>::LPType GaussQueue<ET,false>::true_pop()
{
    if(main_queue.empty())
    {
        ++ (gauss_sieve->number_of_points_sampled);
        ++ (gauss_sieve->number_of_points_constructed);
        return sampler->sample();
    }
    else
    {
        #ifndef USE_REGULAR_QUEUE
        LPType next_point = *(main_queue.top());
        delete main_queue.top();
        #else
        LPType next_point = * (main_queue.front());
        delete main_queue.front();
        #endif // USE_REGULAR_QUEUE
        main_queue.pop();
        return next_point;
    }
}

template<class ET>
typename GaussQueue<ET,false>::LPType* GaussQueue<ET,false>::pop_take_ownership()
{
    if(main_queue.empty())
    {
    ++ (gauss_sieve->number_of_points_sampled);
    ++ (gauss_sieve->number_of_points_constructed);
    LPType *next_point_ptr = new LPType (sampler->sample());
    return next_point_ptr;
    }
    else
    {
        #ifndef USE_REGULAR_QUEUE
        LPType* next_point_ptr = main_queue.top();
        #else
        LPType* next_point_ptr = main_queue.front();
        #endif // USE_REGULAR_QUEUE
        main_queue.pop(); //remove pointer from queue.
        return next_point_ptr;
    }
}


template<class ET>
void GaussQueue<ET,false>::push(LPType const & val)
{
    LPType * new_lp = new LPType (val);
    main_queue.push(new_lp);
}

template<class ET>
void GaussQueue<ET,false>::push(LPType && val)
{
    LPType * new_lp = new LPType (std::move(val) );
    main_queue.push(new_lp);
}

template<class ET>
void GaussQueue<ET,false>::give_ownership(LPType * const valptr)
{
    main_queue.push(valptr);
}

template<class ET>
typename GaussQueue<ET,false>::QueueType::size_type GaussQueue<ET,false>::size() const
{
    return main_queue.size();
}

template<class ET>
bool GaussQueue<ET,false>::empty() const
{
    return main_queue.empty();
}

template<class ET>
GaussQueue<ET,false>::~GaussQueue()
{
//TODO: Delete sampler if owned.
    while(! main_queue.empty() )
    {
        #ifndef USE_REGULAR_QUEUE
        delete main_queue.top();
        #else
        delete main_queue.front();
        #endif // USE_REGULAR_QUEUE
        main_queue.pop();
    }
}


#endif // GAUSS_QUEUE_H
