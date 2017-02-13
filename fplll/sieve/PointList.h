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
#include "assert.h"

//Class for (weakly?) sorted list of lattice points.
//includes thread-safe variant(s). May need experiments which implementation is best. (mutex on whole structure on every write, lock-free,...)
//Note that the application is rather lenient on requirement:
// -- keeping the list sorted is not strictly required and can be worked around. It should be _somewhat_ sorted for efficiency reasons.
//Reading and actually using vectors that are marked-for-deletion does not hurt too much.
//The structure of the algorithm allows for relatively simple garbage collection.


template <class ET>
using PointListSingleThreaded = std::forward_list<LatticePoint <ET> >; //list or forward_list? 
									// Where is this ever used? In the SieveClass for the SingleThreaded implementation we are using list<LPType>

//Note: PointListMultiThreaded owns all its lattice vectors that are reachable by forward iteration.

template <class DT> class ListMultiThreaded;
template <class DT> class ListMTNode;
template <class DT> class MTListIterator;

template <class ET>
using PointListMTNode = ListMTNode<LatticePoint<ET> >;
template <class ET>
using PointListIterator=MTListIterator< LatticePoint<ET> >;
template <class ET>
using PointListMultiThreaded=ListMultiThreaded< LatticePoint <ET> >;

//TODO:ListBin
template <class DT>
//using GarbageBin = std::queue< ListMTNode<DT> * >;
class GarbageBin : public std::queue< ListMTNode<DT> * >
{
    public:
    GarbageBin() = default;
    GarbageBin(GarbageBin const &old) = delete;
    GarbageBin(GarbageBin &&old) = default;
    GarbageBin& operator=(GarbageBin const &old) = delete;
    GarbageBin& operator=(GarbageBin &&old) = default;
    ~GarbageBin()
    {
     empty_trash();
    }
    void empty_trash()
    {
    while(!this->empty())
    {
     delete this->front();
     this->pop();
    }
    };
};

template <class DT>
class ListMultiThreaded{
//friend PointListMTNode<ET>;
public:
    using Node=ListMTNode<DT>;
    using DataType    = DT;
    using DataPointer = DT*;
    using AtomicDataPointer = std::atomic<DT *>;
    using NodePointer = ListMTNode<DT> *;
    using AtomicNodePointer = std::atomic<ListMTNode<DT> *>;
    using Iterator=MTListIterator<DT>;
    //using List = ListMultiThreaded<DT>;
    explicit ListMultiThreaded() : //called when only one thread is running
        mutex_currently_writing(),
        start_sentinel_node (new Node),
        end_sentinel_node   (new Node)
        {
        start_sentinel_node->next_node=end_sentinel_node;
        end_sentinel_node->prev_node=start_sentinel_node;
        start_sentinel_node->nodestatus=static_cast<int>(Node::StatusBit::is_first_node);
        end_sentinel_node->nodestatus  =static_cast<int>(Node::StatusBit::is_last_node);
        };
    ListMultiThreaded(ListMultiThreaded const &old)=delete; //No copying via constructor! (due to mutex)
    ListMultiThreaded(ListMultiThreaded && old)=delete; //No moving (mutex)
    ListMultiThreaded & operator=(ListMultiThreaded const & old) = delete; //No copy assignment. (due to mutex)
    ListMultiThreaded & operator=(ListMultiThreaded && old)=delete; //dito

 //TODO: Constructor from SingleThreaded variant.
  ~ListMultiThreaded() //called when only one (master) thread is running
  {
      Iterator next=begin();
      Iterator cur (start_sentinel_node);
      while(!next.is_end())
      {
        delete cur.p;
        cur=next;
        ++next;
      }
      delete cur.p;
      delete next.p;
      #if 0 //calls ++it on end_sentinel_node, avoid
      for(Iterator it(start_sentinel_node); it.p!=nullptr; delete ( (it++).p) ); //must not initialize with it=begin(), since this skips start sentinel.
      #endif
  }

  Iterator begin() const {return start_sentinel_node->next_node.load(std::memory_order_acquire);}; //returns nullptr on empty list. Note that users never see the start sentinel.
  Iterator end() const {return end_sentinel_node;};
  void unlink(Iterator const &pos, GarbageBin<DT> &gb);
  void insert(Iterator const &pos, DT const &val){DT * const tmp = new DT(val); enlist(pos,tmp);}; //inserts a copy of DT just before pos.
  void enlist(MTListIterator<DT> const &pos, DT * const &valref); //moves *DT just before pos, transfering ownership to the list.

private:
  //marked for deletion.
  std::mutex mutex_currently_writing;
  NodePointer const start_sentinel_node; //node before the start of the list. This is never modified, so no atomic here.
  NodePointer const end_sentinel_node; //node after the end of the list, could probably do with a single sentinel.
};
//
//
//

//restriction: Iterator itself is thread-local.

template<class DT>
class MTListIterator{
    public:
    friend ListMultiThreaded<DT>;
    friend void swap(MTListIterator &A, MTListIterator &B){std::swap(A.p,B.p);};
    typedef std::forward_iterator_tag iterator_category;
    using Node=ListMTNode<DT>;
    using DataType    = DT;
    using DataPointer = DT*;
    //MTListIterator(): p(nullptr){};
    MTListIterator() = delete; //should always init with valid object.
    MTListIterator(Node * const & _p): p(_p){assert(_p!=nullptr);};
    MTListIterator(MTListIterator<DT> const &old): p(old.p){};
    MTListIterator(MTListIterator<DT> && old)=default;
    ~MTListIterator(){};
    MTListIterator<DT>& operator=(MTListIterator<DT> other) {swap(*this,other);return *this;};
    MTListIterator<DT>& operator++() {p=p->next_node.load(memory_order_acquire);assert(p!=nullptr);return *this;}; //prefix version
    MTListIterator<DT> operator++(int) {auto tmp=p; ++(*this);return tmp;}; //postfix version
    //Note: there is no operator--. This is intentional: Assuming the list is only traversed in one direction makes things easier wrt. concurrency.
    DataPointer operator->() const {return p->latpoint;}; //Note weird semantics of -> overload cause latpoint to get dereferenced as well.
    DataType & operator*() const {return *(p->latpoint);};
    operator DataPointer() const {return p->latpoint;}; //converts to pointer to lattice point.
    bool operator==(MTListIterator<DT> const & other) const {return p==other.p;};
    bool operator!=(MTListIterator<DT> const & other) const {return p!=other.p;};
    bool is_end() const {return p->check_for_end_node();};
    bool is_good() const {return !(p->is_marked_for_deletion() );};
    private:
    Node * p; //does not own. Need not be atomic.
};

//owns the lattice point.
template <class DT>
class ListMTNode{
    friend ListMultiThreaded<DT>;
    friend MTListIterator<DT>;
    using DataType    = DT;
    using DataPointer = DT *;
    using AtomicDataPointer = std::atomic<DT *>;
    using NodePointer = ListMTNode<DT> *;
    using AtomicNodePointer = std::atomic<ListMTNode<DT> *>;
    public:
    ListMTNode() : next_node(nullptr),prev_node(nullptr), latpoint(nullptr),nodestatus(0)  {} ;
    ListMTNode(ListMTNode const &old) = delete;
    ListMTNode(ListMTNode &&old) = delete;
    ListMTNode & operator=(ListMTNode const &old) =delete;
    ListMTNode & operator=(ListMTNode &&old) = delete;
    ~ListMTNode(){delete latpoint;}; //destructor
    bool check_for_end_node() const {return nodestatus == static_cast<int>(StatusBit::is_last_node);};
    bool is_marked_for_deletion() const {return nodestatus == static_cast<int>(StatusBit::is_to_be_deleted);};
    bool is_sentinel_node() const {return (nodestatus == static_cast<int>(StatusBit::is_first_node))||(nodestatus == static_cast<int>(StatusBit::is_last_node) );};
    bool is_plain_node() const {return nodestatus == 0;};
    private:
    AtomicNodePointer next_node;
    NodePointer prev_node;
    DataPointer latpoint; //actual data. We may use a pointer here for potential atomicity. This is hidden from the user.
    int nodestatus; //actually a bitfield (not using std:bitset since it only converts safely to ulong or ulonglong)
    public:
    enum class StatusBit
      {
      is_to_be_deleted=1,
      is_first_node=2,
      is_last_node=4
      };

};


template <class DT>
void ListMultiThreaded<DT>::unlink(Iterator const & pos, GarbageBin<DT> &gb)
{
    assert(pos.p!=nullptr);
    {
        assert(!(pos.p->is_sentinel_node()));
        //std::lock_guard<std::mutex> writelock(mutex_currently_writing);
        mutex_currently_writing.lock();
        if(pos.p->is_plain_node()) //otherwise, already deleted.
            {
                pos.p->nodestatus=static_cast<int>(Node::StatusBit::is_to_be_deleted);
                NodePointer nextpos=pos.p->next_node; //.load(memory_order_relaxed); //All writes are within locks anyway.
                nextpos->prev_node=pos.p->prev_node;
                pos.p->prev_node->next_node.store(nextpos,memory_order_relaxed); //relaxed should be fine!!!
                mutex_currently_writing.unlock();
                //Put in garbage bin. //TODO: More clever garbage bin, requires changing structs and global counters.
                gb.push(pos.p);
            }
            else
            {
            mutex_currently_writing.unlock();
            }
    }
    return;
}



template <class DT>
void ListMultiThreaded<DT>::enlist(MTListIterator<DT> const &pos, DT * const &valref)
{
Node* newnode = new Node;
newnode->latpoint = valref;
Node* nextgood =pos.p; //we work directly with the underlying pointer, not using the iterator, since we do not need atomic loads to traverse here.
assert(nextgood!=nullptr);
mutex_currently_writing.lock();
while(nextgood->is_marked_for_deletion()){nextgood=nextgood->next_node;} //.load(memory_order_relaxed);}
Node* preced=nextgood->prev_node;
newnode->next_node=nextgood;
newnode->prev_node=preced;
nextgood->prev_node=newnode;
//until here, no other thread can observe our writes, since we did not publish the pointer.
preced->next_node.store(newnode,memory_order_release);
//this one can be observed (even prior to releasing the lock), and we have to make sure that other non-mutex-protected threads that see this write also see the (non-atomic) writes to newnode.
mutex_currently_writing.unlock();
}
//template <class ET>
//using PointListMultiThreaded= ListMultiThreaded<LatticePoint<ET>>;
//template <class ET> class PointListIterator;

#endif
