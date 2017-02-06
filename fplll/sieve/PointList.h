//
// PointList.h
//

#ifndef lattice_point_list_class_h
#define lattice_point_list_class_h

#include "LatticePoint.h"
#include "sieve_common.h"
#include <mutex>
#include <atomic>
#include <forward_list>
#include <queue>

//Class for (weakly?) sorted list of lattice points.
//includes thread-safe variant(s). May need experiments which implementation is best. (mutex on whole structure on every write, lock-free,...)
//Note that the application is rather lenient on requirement:
// -- keeping the list sorted is not strictly required and can be worked around. It should be _somewhat_ sorted for efficiency reasons.
//Reading and actually using vectors that are marked-for-deletion does not hurt too much.
//The structure of the algorithm allows for relatively simple garbage collection.


template <class ZT>
using PointListSingleThreaded = std::forward_list<LatticePoint <ZT> >; //list or forward_list?

//Note: PointListMultiThreaded owns all its lattice vectors that are reachable by forward iteration.

template <class ZT> class PointListMultiThreaded;
template <class ZT> class PointListMTNode;
template <class ZT> class PointListIterator;
template <class ZT>
using GarbageBin = std::queue< PointListMTNode<ZT> * >;
//template <class ZT> class PointListConstIterator;

template <class ZT>
class PointListMultiThreaded{
//friend PointListMTNode<ZT>;
public:
  explicit PointListMultiThreaded() :
        mutex_currently_writing(),
        start_sentinel_node (new PointListMTNode<ZT>),
        end_sentinel_node   (new PointListMTNode<ZT>)
        {
        start_sentinel_node->next_node=end_sentinel_node;
        end_sentinel_node->prev_node=start_sentinel_node;
        start_sentinel_node->nodestatus=static_cast<int>(PointListMTNode<ZT>::StatusBit::is_first_node);
        end_sentinel_node->nodestatus  =static_cast<int>(PointListMTNode<ZT>::StatusBit::is_last_node);
        };
  PointListMultiThreaded(PointListMultiThreaded const &old)=delete; //No copying via constructor! (due to mutex)
  PointListMultiThreaded(PointListMultiThreaded && old)=delete; //No moving (mutex)
  PointListMultiThreaded & operator=(PointListMultiThreaded const & old) = delete; //No copy assignment. (due to mutex)
  PointListMultiThreaded & operator=(PointListMultiThreaded && old)=delete; //dito
 //TODO: Constructor from SingleThreaded variant.
  ~PointListMultiThreaded()
  {
     for(PointListIterator<ZT> it=PointListIterator<ZT>(start_sentinel_node); it.p!=nullptr; delete ( (it++).p) );
  } //destructor yet missing, leaking memory.

  PointListIterator<ZT> begin(){return start_sentinel_node->next_node;}; //returns nullptr on empty list
  PointListIterator<ZT> end() {return end_sentinel_node;};
  void unlink(PointListIterator<ZT> const &pos, GarbageBin<ZT> &gb);
private:
  //marked for deletion.
  std::mutex mutex_currently_writing;
  PointListMTNode<ZT>* const start_sentinel_node; //node before the start of the list
  PointListMTNode<ZT>* const end_sentinel_node; //node after the end of the list
};



template<class ZT>
class PointListIterator{
    public:
    friend PointListMultiThreaded<ZT>;
    friend void swap(PointListIterator &A, PointListIterator &B){std::swap(A.p,B.p);};
    typedef std::forward_iterator_tag iterator_category;
    PointListIterator(): p(nullptr){};
    PointListIterator(PointListMTNode<ZT>* const & _p): p(_p){};
    PointListIterator(PointListIterator<ZT> const &old): p(old.p){};
    PointListIterator(PointListIterator<ZT> && old)=default;
    ~PointListIterator(){};
    PointListIterator<ZT>& operator=(PointListIterator<ZT> other) {swap(*this,other);return *this;};
    PointListIterator<ZT>& operator++() {p=p->next_node;return *this;}; //prefix version
    PointListIterator<ZT> operator++(int) {auto tmp=p; ++(*this);return tmp;}; //postfix version
    //Note: there is no operator--. This is intentional: Assuming the list is only traversed in one direction makes things easier wrt. concurrency.
    LatticePoint<ZT>* operator->() const {return p->latpoint;}; //Note weird semantics of -> overload cause latpoint to get dereferenced as well.
    LatticePoint<ZT>& operator*() const {return *(p->latpoint);};
    operator LatticePoint<ZT>*() const {return p->latpoint;}; //converts to pointer to lattice point.
    bool operator==(PointListIterator<ZT> const & other) const {return p==other.p;};
    bool operator!=(PointListIterator<ZT> const & other) const {return p!=other.p;};
    bool is_end() const {return p->check_for_end_node();};
    private:
    PointListMTNode<ZT> * p; //does not own.
};

//owns the lattice point.
template <class ZT>
class PointListMTNode{
    friend PointListMultiThreaded<ZT>;
    friend PointListIterator<ZT>;
    public:
    PointListMTNode() : next_node(nullptr),prev_node(nullptr), latpoint(nullptr),nodestatus(0)  {} ;
    PointListMTNode(PointListMTNode const &old) = delete;
    PointListMTNode(PointListMTNode &&old) = delete;
    PointListMTNode & operator=(PointListMTNode const &old) =delete;
    PointListMTNode & operator=(PointListMTNode &&old) = delete;
    ~PointListMTNode(){delete latpoint;}; //destructor
    bool check_for_end_node() const {return nodestatus == static_cast<int>(StatusBit::is_last_node);};
    bool is_marked_for_deletion() const {return nodestatus == static_cast<int>(StatusBit::is_to_be_deleted);};
    bool is_sentinel_node() const {return (nodestatus == static_cast<int>(StatusBit::is_first_node))||(nodestatus == static_cast<int>(StatusBit::is_last_node) );};
    bool is_plain_node() const {return nodestatus == 0;};
    private:
    PointListMTNode *next_node;
    PointListMTNode *prev_node;
    LatticePoint<ZT> *latpoint; //actual data. We may use a pointer rather than an actual latpoint here to allow atomic replacements. This is hidden from the user.
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
