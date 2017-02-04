//
// PointList.h
//

#ifndef lattice_point_list_class_h
#define lattice_point_list_class_h

#include "LatticePoint.h"
#include "sieve_common.h"
#include <mutex>

//Class for (weakly?) sorted list of lattice points.
//includes thread-safe variant(s). May need experiments which implementation is best. (mutex on whole structure on every write, lock-free,...)
//Note that the application is rather lenient on requirement:
// -- keeping the list sorted is not strictly required and can be worked around. It should be _somewhat_ sorted for efficiency reasons.
//Reading and actually using vectors that are marked-for-deletion does not hurt too much.
//The structure of the algorithm allows for relatively simple garbage collection.


//typedef std::list<typename LatticePoint> PointListSingleThreaded; //technically, forward_list would be good enough.

template <class ZT>
using PointListSingleThreaded = std::list<LatticePoint <ZT> >;

//Note: PointListMultiThreaded owns its lattice vectors.

template <class ZT> class PointListMultiThreaded;
template <class ZT> class PointListMTNode;
template <class ZT> class PointListIterator;
//template <class ZT> class PointListConstIterator;

template <class ZT>
class PointListMultiThreaded{
friend PointListMTNode<ZT>;
public:
  explicit PointListMultiThreaded() :
        start_sentinel_node (new PointListMTNode<ZT>),
        end_sentinel_node   (new PointListMTNode<ZT>),
        mutex_currently_writing()
        {
        start_sentinel_node->next_node=end_sentinel_node;
        end_sentinel_node->prev_node=start_sentinel_node;
        start_sentinel_node->nodestatus=PointListMTNode<ZT>::StatusBit::is_first_node;
        end_sentinel_node->nodestatus=PointListMTNode<ZT>::StatusBit::is_last_node;
        };
  PointListMultiThreaded(PointListMultiThreaded const &old)=delete; //No copying via constructor!
  PointListMultiThreaded(PointListMultiThreaded && old)=default; //moving is OK.
  PointListMultiThreaded & operator=(PointListMultiThreaded const & old) = delete; //No copy assignment.
  PointListMultiThreaded & operator=(PointListMultiThreaded && old)=default; //Move assignment OK.
 //TODO: Constructor from SingleThreaded variant.
  ~PointListMultiThreaded(); //destructor
  PointListIterator<ZT> begin(){return start_sentinel_node->next_node;}; //returns nullptr on empty list
  PointListIterator<ZT> end() {return end_sentinel_node;};
private:
  //marked for deletion.
  std::mutex mutex_currently_writing;
  PointListMTNode<ZT>* const start_sentinel_node; //node before the start of the list
  PointListMTNode<ZT>* const end_sentinel_node; //node after the end of the list
};

template<class ZT>
class PointListIterator{
    public:
    friend void swap(PointListIterator &A, PointListIterator &B){std::swap(A.p,B.p);};
    typedef std::forward_iterator_tag iterator_category;
    PointListIterator(): p(nullptr){};
    PointListIterator(PointListMTNode<ZT>* const & _p): p(_p){};
    PointListIterator(PointListIterator<ZT> const &old): p(old.p){};
    PointListIterator(PointListIterator<ZT> && old)=default;
    ~PointListIterator(){};
    PointListIterator<ZT>& operator=(PointListIterator<ZT> other) {swap(*this,other);return *this;};
    PointListIterator<ZT>& operator++() {p=p->next_node;return *this;}; //prefix version, will give a non-deferen
    PointListIterator<ZT> operator++(int) {auto tmp=p; ++(*this);return tmp;}; //postfix version, segfaults when used on last element.
    //Note: there is no operator--. This is intentional: Assuming the list is only traversed in one direction makes things easier wrt. concurrency.
    LatticePoint<ZT>* operator->() const {return p;}; //Note weird semantics of -> overload.
    LatticePoint<ZT>& operator*() const {return *p;};
    bool operator==(PointListIterator<ZT> const & other) const {return p==other.p;};
    bool operator!=(PointListIterator<ZT> const & other) const {return p!=other.p;};
    private:
    PointListMTNode<ZT> *p;
};

template <class ZT>
class PointListMTNode{
    public:
    PointListMTNode() : next_node(nullptr),prev_node(nullptr), latpoint(nullptr),nodestatus(0)  {} ;
    PointListMTNode(PointListMTNode const &old);
    PointListMTNode & operator=(PointListMTNode const &old);
    ~PointListMTNode(); //destructor
    private:
    PointListMTNode *next_node;
    PointListMTNode *prev_node;
    LatticePoint<ZT> *latpoint; //actual data
    int nodestatus; //actually a bitfield (not using std:bitset since it only converts safely to ulong or ulonglong)
    public:
    enum class StatusBit
      {
      is_to_be_deleted=1,
      is_first_node=2,
      is_last_node=4
      };

};

#endif
