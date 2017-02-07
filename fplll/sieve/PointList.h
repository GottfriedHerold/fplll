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

template <class DT> class ListMultiThreaded;
template <class DT> class ListMTNode;
template <class DT> class MTListIterator;

//template <class ZT> class PointListMultiThreaded;
//template <class ZT> class PointListMTNode;
template <class ZT>
using PointListMTNode = ListMTNode<LatticePoint<ZT> >;
template <class ZT>
using PointListIterator=MTListIterator< LatticePoint<ZT> >;
template <class ZT>
using PointListMultiThreaded=ListMultiThreaded< LatticePoint <ZT> >;

//TODO:ListBin
template <class DT>
using GarbageBin = std::queue< ListMTNode<DT> * >;

template <class DT>
class ListMultiThreaded{
//friend PointListMTNode<ZT>;
public:
    using Node=ListMTNode<DT>;
    using DataType    = DT;
    using DataPointer = DT*;
    using Iterator=MTListIterator<DT>;
    //using List = ListMultiThreaded<DT>;
    explicit ListMultiThreaded() :
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
  ~ListMultiThreaded()
  {
      for(Iterator it(start_sentinel_node); it.p!=nullptr; delete ( (it++).p) ); //must not initialize with it=begin(), since this skips start sentinel.
  } //destructor yet missing, leaking memory.

  Iterator begin(){return start_sentinel_node->next_node;}; //returns nullptr on empty list. Note that users never see the start sentinel.
  Iterator end() {return end_sentinel_node;};
  void unlink(Iterator const &pos, GarbageBin<DT> &gb);
  void insert(Iterator const &pos, DT const &val){DT * const tmp = new DT(val); enlist(pos,tmp);}; //inserts a copy of DT just before pos.
  void enlist(MTListIterator<DT> const &pos, DT * const &valref); //moves *DT just before pos, transfering ownership to the list.

private:
  //marked for deletion.
  std::mutex mutex_currently_writing;
  Node* const start_sentinel_node; //node before the start of the list
  Node* const end_sentinel_node; //node after the end of the list
};
//
//
//
template<class DT>
class MTListIterator{
    public:
    friend ListMultiThreaded<DT>;
    friend void swap(MTListIterator &A, MTListIterator &B){std::swap(A.p,B.p);};
    typedef std::forward_iterator_tag iterator_category;
    using Node=ListMTNode<DT>;
    using DataType    = DT;
    using DataPointer = DT*;
    MTListIterator(): p(nullptr){};
    MTListIterator(Node * const & _p): p(_p){};
    MTListIterator(MTListIterator<DT> const &old): p(old.p){};
    MTListIterator(MTListIterator<DT> && old)=default;
    ~MTListIterator(){};
    MTListIterator<DT>& operator=(MTListIterator<DT> other) {swap(*this,other);return *this;};
    MTListIterator<DT>& operator++() {p=p->next_node;return *this;}; //prefix version
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
    Node * p; //does not own.
};

//owns the lattice point.
template <class DT>
class ListMTNode{
    friend ListMultiThreaded<DT>;
    friend MTListIterator<DT>;
    using DataType    = DT;
    using DataPointer = DT *;
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
    ListMTNode *next_node;
    ListMTNode *prev_node;
    DataPointer latpoint; //actual data. We may use a pointer rather than an actual latpoint here to allow atomic replacements. This is hidden from the user.
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
    {
        assert(!(pos.p->is_sentinel_node()));
        //std::lock_guard<std::mutex> writelock(mutex_currently_writing);
        mutex_currently_writing.lock();
        if(pos.p->is_plain_node()) //otherwise, already deleted.
            {
                pos.p->nodestatus=static_cast<int>(Node::StatusBit::is_to_be_deleted);
                pos.p->next_node->prev_node=pos.p->prev_node;
                pos.p->prev_node->next_node=pos.p->next_node;
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
Iterator next_good(pos);
for (;;++next_good)
{
if (!next_good.is_good())
    {
    continue;
    }
mutex_currently_writing.lock();
if (next_good.is_good())
    {
    break;
    }
mutex_currently_writing.unlock();
}
newnode->next_node=next_good.p;
newnode->prev_node=next_good.p->prev_node;
newnode->next_node->prev_node=newnode;
newnode->prev_node->next_node=newnode;
//does not work -- need atomics to prevent compiler from reorderings.
mutex_currently_writing.unlock();
}

//template <class ZT>
//using PointListMultiThreaded= ListMultiThreaded<LatticePoint<ZT>>;
//template <class ZT> class PointListIterator;




#endif
