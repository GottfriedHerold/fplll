#ifndef ARRAY_LIST_H
#define ARRAY_LIST_H

#include "DefaultIncludes.h"

namespace GaussSieve
{

/**
  ArrayList is an list of arrays that offer a compromise between arrays (consecutive in memory),
  which are optimized for fast traversal, and lists, wich are optimized for arbitrary insertion.

  ArrayList<T,blocksize> is realized as circularly doubly-linked list of blocks with a special
  sentinel block (that acts as both before-begin and after-the-end position of the data structure).
  Each block is responsible for up to blocksize many elements of type T.
  The number of elements actually used is stored in used_size, which is guaranteed to be >0 iff and
  only if the block is not the sentinel block.
  The actual data are stored in a (dynamically allocated, using placement-new in preallocated
  storage) c-style array whose address is stored in each block.
  (Note: iterators maintain a copy of that pointer to avoid resolving an extra level of indirection
  during traversal)
  During insertion / deletions, we create or merge blocks as needed.
*/

namespace ArrayListDetails
{
// Our doubly-linked list nodes (i.e. blocks of the array list) consist of 2 parts:
// The linked-list (meta)data structure and the payload data (a cstyle-array/pointer in our case)
// We use a generic head class for the meta-data; to create the node class, we inherit from it.
// The advantage is that this setup is that the meta-data does not depend on template paramters.
// The next / prev pointer only point to the base class.
// As a consequence, we can in principle have linked lists whose payload data type differs among
// nodes; in particular, we may use a different payload for the sentinel nodes.
// This way, the sentinel node actually IS the ArrayList and its payload consists of additional
// data like total size etc.
// (this is not so relevant if the payload is just a pointer, but I want to be able to experiment
// with keeping the data directly inside the nodes, in which case it makes a difference if
// totalsize/blocksize is small. In particular, the efficiency should actually match that of a
// plain c-style array in the limiting case.)
class LinkedListNodeMeta;

class LinkedListNodeMeta
{
  LinkedListNodeMeta* next;  // pointer to next block.
  LinkedListNodeMeta* prev;  // pointer to previous block.
  unsigned int used_size;     // number of actual data used in the block. 0 for the sentinel.
                             // (so also act as an indicator whether this is the sentinel block.)
  explicit LinkedListNodeMeta(LinkedListNodeMeta * const new_next = nullptr,
                              LinkedListNodeMeta * const new_prev = nullptr,
                              unsigned int const new_used_size = 0)
  : next(new_next), prev(new_prev), used_size(new_used_size)
  {
  }
};

template<class T, unsigned int blocksize>
class ArrayListBlock : LinkedListNodeMeta
{
  static_assert(blocksize > 0, "");
  T *data;  // C-array
};

}  // end namespace ArrayListDetails

template<class T, unsigned int blocksize>
class ArrayListConstIterator
{
  unsigned int index; // index inside block
  T *dataptr;  // dataptr (nullptr for end-iterator)
  ArrayListDetails::LinkedListNodeMeta* nodeptr; // pointer to relevant node
  static_assert(blocksize > 0, "");
};

template<class T, unsigned int blocksize>
class ArrayList: ArrayListDetails::LinkedListNodeMeta
{
  static_assert(blocksize > 0, "");
  // creates an empty array list

public:

  using value_type = T;
  using size_type = std::size_t;
  using const_iterator = ArrayListConstIterator<T, blocksize>;
  using reference = T&;
  using const_reference = T const &;

  explicit ArrayList() noexcept : LinkedListNodeMeta(this, this, 0), total_size(0) {}
  ArrayList(ArrayList const &) = delete;
  ArrayList(ArrayList &&) = default;

  reference front() = delete;
  const_reference front() const = delete;
  reference back()  = delete;
  const_reference back() const = delete;
  const_iterator cbegin() const noexcept = delete;
  const_iterator cend()   const noexcept = delete;
  NODISCARD bool empty() const noexcept = delete;
  size_type size() const noexcept { return total_size; }
  void clear() noexcept = delete;
  // insert and insert_after
  // emplace ???
  // erase
  // push_back
  // emplace_back
  // pop_back : replace by true_pop_back
  // sort


private:
  size_type total_size;
  using Base = ArrayListDetails::LinkedListNodeMeta;

  Base* insert_block_between(Base &before, Base &after)
  {
#ifdef DEBUG_SIEVE_ARRAYLIST
    assert(before.next == after);
    assert(after.prev  == before);
#endif

  }
};



//namespace ArrayListDetails
//{
//
///**
//  ArrayListBlock<T,blocksize> is one single block of an array list.
//  We use a fixed-size array of blocksize T's for that.
//  Note that typically, not all entries of a given block are actually used.
//  Rather, the used blocks entries are given by
//  block[usedsize-1],...,block[0].
//  Note: For optimization reasons, the index inside the block is counting down as we traverse forward
//        The reason is that we optimize for forward traversal and this avoids comparing against (and
//        needing to re-read) usedsize to determine whether we are at the end of a block.
//  Note that we will never have empty blocks.
//*/
//template<class T, unsigned int blocksize>
//struct ArrayListBlock
//{
//  std::array<T,blocksize> block;
//  unsigned int usedsize;
//
//  ArrayListBlock() : block(), usedsize(0) {}  // Note: This creates an empty block, which needs to
//                                              // be filled with at least one element.
//  ArrayListBlock(T const &value): block(), usedsize(1)
//  {
//    block[0] = value;
//  }
//  ArrayListBlock(T &&value): block(), usedsize(1)
//  {
//    block[0] = std::move(value);
//  }
//
//  ArrayListBlock(ArrayListBlock const &) = default;
//  ArrayListBlock(ArrayListBlock &&) = default;
//};
//
//}  // end namespace ArrayListDetails
//
//template<class T, unsigned int blocksize> class ArrayList;
//template<class T, unsigned int blocksize> class ArrayListConstIterator;
//
//template<class T, unsigned int blocksize>
//class ArrayList
//{
//  static_assert(blocksize > 0, "");
//  static_assert(std::is_default_constructible<T>::value, "");
//
//private:
//  using Block = ArrayListDetails::ArrayListBlock<T,blocksize>;
//  std::list<Block> blocks;
//  std::size_t total_size; // total amount of actually used entries. Note that this disregards a
//                          // "dummy" before-begin entry that we always keep.
//public:
//  using const_iterator = ArrayListConstIterator<T,blocksize>;
//  using size_type = std::size_t;
//  using value_type = T;
//  using block_iterator = typename decltype(blocks)::iterator;
//  friend const_iterator;
//
//  ArrayList() : blocks(), total_size(0)
//  {
//    // creates a dummy before-begin entry. It purpose is mostly to eliminate some special cases.
//    blocks.emplace_back();
//    blocks.front().usedsize = 1;
//  }
//  [[deprecated("Copying array lists is probably wrong")]]ArrayList(ArrayList const &) = default;
//  ArrayList(ArrayList &&) = default;  // Note: This invalidates iterators.
//
//  size_type size() const noexcept { return total_size; }
//  NODISCARD bool empty() const noexcept { return total_size == 0; }
//  void clear() noexcept
//  {
//    blocks.clear();
//    total_size = 0;
//    // Create a dummy entry at the beginning (which does not contribute to total_size)
//    blocks.emplace_back();
//    blocks.front().usedsize = 1;
//  }
//
////  block_iterator insert_empty_block_before(block_iterator const &block_it)
////  {
////    return blocks.emplace(block_it);
////  }
//
//private:
//  // internal management function:
//  // The argument it is an iterator to an overfull block.
//  // The function will move parts of the overfull block into adjancent blocks
//  // (possibly creating a new adjacent block for that purpose)
//  // This invalidates iterators. The return value points to the same object as the input.
//  const_iterator split_blocks(const_iterator const &it)
//  {
//    block_iterator const current_block_it = it.list_it;
//    assert(current_block_it->usedsize == blocksize); // for now
//    assert(false);
//  }
//
//public:
//  const_iterator cbegin() const noexcept
//  {
//    const_iterator ret{blocks.cbegin(),blocks.cend()};
//    ++ret;
//    return ret;
//  }
//
//  const_iterator cbefore_begin() const noexcept
//  {
//    return const_iterator{blocks.cbegin(), blocks.cend()};
//  }
//
//  const_iterator cend() const noexcept
//  {
//    return const_iterator{blocks.cend(), 0, blocks.cend()};
//  }
//
//
//
////  void push_back(T const &value)
////  {
////
////  }
//
//
//
//};
//
//template<class T, unsigned int blocksize>
//class ArrayListConstIterator
//{
//  friend ArrayList<T, blocksize>;
//private:
//  using Block = ArrayListDetails::ArrayListBlock<T,blocksize>;
//  using Array = decltype (Block::block);
//  using ArrayIt  = typename Array::iterator;
//  using ListIt   = typename std::list<Block>::iterator;
//  using MySelf   = ArrayListConstIterator<T,blocksize>;
//
//  ListIt list_it;
//  unsigned int index;
//  ListIt end_it;  // could make it a const_iterator
////  unsigned int current_block_size;
//
////  void update_max_index() { current_block_size = list_it->used_size; }
//
//  // Default constructors and destructors
//public:
//
//  bool operator==(MySelf const &other) const { return (list_it==other.list_it) && (index==other.index); }
//  bool operator!=(MySelf const &other) const { return (list_it!=other.list_it) || (index!=other.index); }
//
//  MySelf operator=(MySelf const &) = default;
//  MySelf operator=(MySelf &&) = default;
//
//  ArrayListConstIterator(ListIt const &new_list_it, unsigned int const &new_index, ListIt const &new_end_it)
//      : list_it(new_list_it), index(new_index), end_it(new_end_it) {}
//
//  // points to the beginning of the block indicated by new_list_it
//  ArrayListConstIterator(ListIt const &new_list_it, ListIt const &new_end_it)
//      : list_it(new_list_it), index(new_list_it->used_size-1), end_it(new_end_it)
//  {
//    assert(new_list_it->used_size>0);
//  }
//
//  MySelf &operator++()  // prefix only for now.
//  {
//    if(index>0)
//    {
//      --index;
//    }
//    else
//    {
//      ++list_it;
//      if (list_it == end_it)
//      {
//        index = 0;
//      }
//      else
//      {
//        index = list_it->used_size - 1;
//      }
//    }
//    return *this;
//  }
//
//  T const & operator*() const { return list_it->block[index]; }
//  T const * operator->() const { return list_it->block.data()+index; }
//
//};
//


}  // end namespace GaussSieve


#endif
