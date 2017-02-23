//
// PointList.h
//

#ifndef lattice_point_list_class_h
#define lattice_point_list_class_h


//ET : Entry tpe : type of entries of vectors we are dealing with. Usually an integral type like ET = Z_NR<long> or Z_NR<mpz_t>
//DT : (fundamental) Data type : entries in our custom containers, i.e. lattice points, usually  DT = LatticePoint<ET>


//forward declarations.

template <class ET, bool MT, int n_fixed = -1>
class GaussList; //holds approximations

template <class ET, bool MT, int n_fixed = -1>
class GaussIterator;
//template<class ET>
//using PointListSingleThreaded = GaussList<ET, false>;

//template <class ET>
//using PointListSingleThreaded = std::forward_list<LatticePoint <ET> >; //list or forward_list?
//
//Note: PointListMultiThreaded owns all its lattice vectors that are reachable by forward iteration.

template <class DT> class ListMultiThreaded;
template <class DT> class ListMTNode;
template <class DT> class MTListIterator;
//
//template <class ET>
//using PointListMTNode = ListMTNode<LatticePoint<ET> >;
//template <class ET>
//using PointListIterator=MTListIterator< LatticePoint<ET> >;
//template <class ET>
//using PointListMultiThreaded=ListMultiThreaded< LatticePoint <ET> >;
//
template <class DT> class GarbageBin;


#include "LatticePoint.h"
#include "sieve_common.h"
#include <mutex>
#include <atomic>
#include <forward_list>
#include <queue>
#include "assert.h"
#include "LatticePoint2.h"

//Class for (weakly?) sorted list of lattice points.
//includes thread-safe variant(s). May need experiments which implementation is best. (mutex on whole structure on every write, lock-free,...)
//Note that the application is rather lenient on requirement:
// -- keeping the list sorted is not strictly required and can be worked around. It should be _somewhat_ sorted for efficiency reasons.
//Reading and actually using vectors that are marked-for-deletion does not hurt too much.
//The structure of the algorithm allows for relatively simple garbage collection.
//Note that we do NOT assume a sequentially constistent memory model here. Relaxing from that should gain a bit of performance:
//Be aware that reading time from the lattice point list even asymptotically leading order.

template <class ET>
class GaussList<ET, false, -1>
{
public:
    using EntryType= ET;
    using DataType = ApproxLatticePoint<ET,false,-1>;
    using DataPointer=DataType *;
    using UnderlyingContainer = std::list<DataType>;
    //using Iterator = typename UnderlyingContainer::const_iterator;
    using Iterator = GaussIterator<ET,false,-1>;
    using DetailType = typename DataType::DetailType;
    //using NonConstIterator = typename UnderlyingContainer::iterator;
    explicit GaussList() = default;
    GaussList(GaussList const & old) = delete;
    GaussList(GaussList && old) = delete;
    GaussList & operator= (GaussList const &other) = delete;
    GaussList & operator= (GaussList &&other) = delete;
    ~GaussList() = default; //just deletes the underlying container. Fine as long as we don't store pointers.
    Iterator cbegin()                                                   {return actual_list.begin() ;};
    Iterator cend()                                                     {return actual_list.end();};
    //Iterator cbefore_begin() const                                    {return actual_list.cbefore_begin();};
    Iterator insert_before(Iterator pos, DetailType const & val)        {return actual_list.emplace(pos, val);};
    //Iterator insert_before(Iterator pos, DetailType && val)           {return actual_list.insert(pos,std::move(val));};
    Iterator erase(Iterator pos)                                        {return actual_list.erase(pos);};
    //void unlink(Iterator pos)                                           {actual_list.erase(pos);}; //only for single-threaded
    void sort()                                                         {actual_list.sort();};  //only for single-threaded (for now)

private:
    UnderlyingContainer actual_list;
};

template <class ET>
class GaussIterator<ET,false,-1> : public std::list<ApproxLatticePoint<ET,false, -1 > >::iterator
{
    public:
    LatticePoint<ET> * access_details() { assert(  (*this)->get_details_ptr_rw()!=nullptr );  return (*this)-> get_details_ptr_rw() ;};
    GaussIterator(typename std::list< ApproxLatticePoint<ET,false, -1 > >::iterator other) : std::list<ApproxLatticePoint<ET,false,-1> >::iterator(other) {};
};



//thread-safe lock-free list.
//Ownership: All points not on some garbage queue are owned by the list.
//Important limitations:
//YOU MAY NOT EDIT ELEMENTS VIA ITERATORS. Make a copy of (*it) and modify the copy. Delete the old element from the list and put a new one into it.
//The destructor is not thread-safe. Make sure only one master thread is creating / deleting the lists.
//When deleting, no other thread must access the list any longer. Garbage cleaning is independent from deleting the list.
//Iterators should be thread-local. Do not pass iterators between threads.
//Iterators are forward iterators (i.e. can only be increased)
//All iterators are const_iterators (i.e. trying to edit *it will give a compile-time error).


template <class DT>
class ListMultiThreaded
{
//friend PointListMTNode<ET>;
public:
    using Node=ListMTNode<DT>;
    using DataType    = DT;
    using DataPointer = DT*;
    using AtomicDataPointer = std::atomic<DT *>;
    using NodePointer = ListMTNode<DT> *;
    using AtomicNodePointer = std::atomic<ListMTNode<DT> *>;
    using Iterator=MTListIterator<DT>;
    explicit ListMultiThreaded(); //called when only one thread is running. Creates empty list.
    ListMultiThreaded(ListMultiThreaded const &old)=delete; //No copying via constructor! (due to mutex)
    ListMultiThreaded(ListMultiThreaded && old)=delete; //No moving (mutex)
    ListMultiThreaded & operator=(ListMultiThreaded const & old) = delete; //No copy assignment. (due to mutex)
    ListMultiThreaded & operator=(ListMultiThreaded && old)=delete; //dito
    //TODO: Constructor from SingleThreaded variant and vice versa?
    ~ListMultiThreaded();
    Iterator cbegin() const; //returns iterator to first element. On empty list, returns valid past-the-end iterator.
    Iterator cend() const;   //returns past-the-end iterator. Must not be dereferenced. (same behaviour as stl::list)
    Iterator cbefore_begin() const; //returns iterator before the first element. May not be dereferenced. Included for completeness.
    void unlink(Iterator const &pos, GarbageBin<DT> &gb); //marks the element pointed at by pos as deleted. Such elements will (eventually)
                                                          //become unreachable by traversing the list. Any Iterators to *pos (including pos)
                                                          //remain valid. *pos is put onto gb, whose job is eventually freeing memory.
                                                          //TODO : Allow forwarding additional arguments to gb.
        //  For each of the following insertion routines, the
        //  return value is an iterator to the newly inserted element.
    Iterator insert_before(Iterator const &pos, DT const &val);     //inserts a copy of val before pos. The iterator pos remains valid.
                                                                    //NOTE: If *pos is marked-for-deletion, val will be inserted before
                                                                    //the next non-marked position pos' that is reachable by increasing pos.
                                                                    //Even in this case, incrementing pos will eventually reach pos' without seeing the newly inserted value.
    Iterator enlist_before(Iterator const &pos, DT * const valref); //moves *valref just before pos (as above), transfering ownership to the list. Avoids copying of *valref.
//same as above, but inserts after the corresponding element. Return values is an iterator to the newly inserted element.
//In case *pos is marked-for-deletion, inserts after the next non-marked position pos' that is reachable by increasing pos.
//Due to these reasons, it is not guaranteed that pos+1 == retval.
//Furthermore, insert_before(pos+1,val) may differ from insert_after(pos,val).
    Iterator insert_after(Iterator  const &pos, DT const &val) = delete;     //not implemented yet.
    Iterator enlist_after(Iterator  const &pos, DT * const valref) = delete; //not implemented yet.

private:
    std::mutex mutex_currently_writing;
    NodePointer const start_sentinel_node; //node before the start of the list. This is never modified, so no atomic here.
    NodePointer const end_sentinel_node; //node after the end of the list, could probably do with a single sentinel.
};

//restriction: Iterator itself is thread-local.

template <class DT> //need to redo. //TODO: Thou shalt not inherit from STL containers (no virtual destructors).
class GarbageBin : public std::queue< ListMTNode<DT> * >
{
public:
    GarbageBin() = default;
    GarbageBin(GarbageBin const &old) = delete;
    GarbageBin(GarbageBin &&old) = default;
    GarbageBin& operator=(GarbageBin const &old) = delete;
    GarbageBin& operator=(GarbageBin &&old) = default;
    ~GarbageBin()   {empty_trash();}
    void empty_trash();
};


template<class DT>
class MTListIterator
{
public:
    friend ListMultiThreaded<DT>;
    friend void swap(MTListIterator &A, MTListIterator &B)      {std::swap(A.p,B.p);};
    using iterator_category = std::forward_iterator_tag;
    using Node=ListMTNode<DT>;
    using DataType    = DT;
    using DataPointer = DT*;
    using CDataPointer= DT const *;
    MTListIterator() = delete; //should always init with valid object.
    MTListIterator(Node * const & _p): p(_p)                    {assert(_p!=nullptr);};
    MTListIterator(MTListIterator<DT> const &old): p(old.p)     {};
    MTListIterator(MTListIterator<DT> && old)=default;
    ~MTListIterator() {};
    MTListIterator<DT>& operator=(MTListIterator<DT> const &other) = default;
    MTListIterator<DT>& operator=(MTListIterator<DT> &&other) = default;
    MTListIterator<DT>& operator++(); //prefix version
    MTListIterator<DT>  operator++(int); //postfix version
    DataPointer operator->() const                              {return p->latpoint;}; //Note weird semantics of -> overload cause latpoint (which is a pointer) to get dereferenced as well.
    DataType & operator*()   const                              {return *(p->latpoint);};
    operator DataPointer() const                                {return p->latpoint;}; //converts from Iterator to pointer to lattice point.
    bool operator==(MTListIterator<DT> const & other) const     {return p==other.p;};  //test for equality of iterators.
    bool operator!=(MTListIterator<DT> const & other) const     {return p!=other.p;};  //test for inequality of iterators. Note: Even if it1 != it2, it may still be the case that *it1 == *it2.
    bool is_end() const                                         {return p->check_for_end_node();};
    bool is_good() const                                        {return !(p->is_marked_for_deletion() ); //Not: Note reliable if caller did not lock mutex. May be made private.
    };
private:
    Node * p; //does not own. Need not be atomic. Is NOT a const - pointer!
};
//TODO: define iterator_traits specialisation (or derive from std::iterator)


//owns the lattice point. For internal use only.
template <class DT>
class ListMTNode
{
    friend ListMultiThreaded<DT>;
    friend MTListIterator<DT>;
    using DataType    = DT;
    using DataPointer = DT*;
    using AtomicDataPointer = std::atomic<DT *>;
    using NodePointer = ListMTNode<DT> *;
    using AtomicNodePointer = std::atomic<ListMTNode<DT> *>;
public:
    ListMTNode() : next_node(nullptr),prev_node(nullptr), latpoint(nullptr),nodestatus(0)  {} ;
    ListMTNode(ListMTNode const &old) = delete;
    ListMTNode(ListMTNode &&old) = delete;
    ListMTNode & operator=(ListMTNode const &old) =delete;
    ListMTNode & operator=(ListMTNode &&old) = delete;
    ~ListMTNode()                                               {delete latpoint;}; //destructor
    bool check_for_end_node() const                             {return nodestatus == static_cast<int>(StatusBit::is_last_node);};
    bool is_marked_for_deletion() const                         {return nodestatus == static_cast<int>(StatusBit::is_to_be_deleted);};
    bool is_sentinel_node() const                               {return (nodestatus == static_cast<int>(StatusBit::is_first_node))||(nodestatus == static_cast<int>(StatusBit::is_last_node) );};
    bool is_plain_node() const                                  {return nodestatus == 0;};
private:
    AtomicNodePointer next_node;
    NodePointer       prev_node;
    DataPointer       latpoint; //actual data. We use a pointer here for potential atomicity. This is hidden from the user.
    int nodestatus; //We use int rather than an enum-type to be on the safe side regarding that atomic-operations be possible.
public:
    enum class StatusBit
    {
        is_to_be_deleted=1,
        is_first_node=2,
        is_last_node=4
    }; //meaning of values for nodestatus. Comparison requires typecast.
};

//template <class ET>
//using PointListMultiThreaded= ListMultiThreaded<LatticePoint<ET>>;
//template <class ET> class PointListIterator;

#include "PointList.cpp"

#endif
