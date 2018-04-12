/**
  Implementation file for ArrayList.h, included from it.
*/

#ifndef ARRAY_LIST_IMPL_H
#define ARRAY_LIST_IMPL_H

namespace GaussSieve
{

namespace ArrayListDetails
{

/**
  methods for Meta and Block
*/

// constructor of block
template<class T, unsigned int blocksize>
template<class... Args>
inline constexpr ArrayListBlock<T,blocksize>::ArrayListBlock(Args&&... args)
  : LinkedListNodeMeta(std::forward<Args>(args)...)
{
  memory_buf = std::malloc(sizeof(T)*blocksize); // Note: C++11 does not have aligned versions.
}

// destructor of block
template<class T, unsigned int blocksize>
inline ArrayListBlock<T, blocksize>::~ArrayListBlock()
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(used_size > 0);
  assert(used_size <= blocksize);
#endif
  // Note: This will destruct objects constructed in-place, but NOT in the order of construction
  // (which we do not know about, anyway)
  for(unsigned int i = 0; i < used_size; ++i)
  {
    reinterpret_cast<T*>(memory_buf)[i].~T();
  }
  std::free(memory_buf);
}

} // end namespace ArrayListDetails (inside namespace GaussSieve)

/**
  Methods for the actual ArrayList
*/

/**
  Methods for the iterator
*/

template<class T, unsigned int blocksize>
inline ArrayListConstIterator<T, blocksize>::ArrayListConstIterator(Meta* const new_nodeptr) noexcept
: nodeptr(new_nodeptr)
{
  if (nodeptr->used_size <= blocksize)
  {
#ifdef DEBUG_SIEVE_ARRAYLIST
    assert(nodeptr->used_size > 0);
#endif
    dataptr = reinterpret_cast<T*> (static_cast<Block*>(nodeptr)->memory_buf);
    index   = nodeptr->used_size - 1;
  }
  else // past-the end or before-begin iterator
  {
#ifdef DEBUG_SIEVE_ARRAYLIST
    assert(nodeptr->used_size == 2*blocksize + 1);
#endif
    dataptr = nullptr;
    index   = 0;
  }
}

template<class T, unsigned int blocksize>
inline ArrayListConstIterator<T, blocksize>::ArrayListConstIterator(Meta* const new_nodeptr, unsigned int const new_index) noexcept
: index(new_index),
  dataptr(reinterpret_cast<T*>(static_cast< Block* >(new_nodeptr)->memory_buf)),
  nodeptr(new_nodeptr)
{ }

template<class T, unsigned int blocksize>
constexpr inline ArrayListConstIterator<T, blocksize>::ArrayListConstIterator(Meta* const new_nodeptr, unsigned int const new_index, T * const new_dataptr) noexcept
: index(new_index), dataptr(new_dataptr), nodeptr(new_nodeptr)
{ }

template<class T, unsigned int blocksize>
inline ArrayListConstIterator<T,blocksize> &ArrayListConstIterator<T, blocksize>::operator++()
{
  if(index > 0) // not at end of block.
  {
    --index; // recall that index counts downwards
    return *this;
  }
  else // jump to next block
  {
    *this = ArrayListConstIterator{nodeptr->next};
    return *this;
  }
}

template<class T, unsigned int blocksize>
inline ArrayListConstIterator<T, blocksize> ArrayListConstIterator<T, blocksize>::operator++(int)
{
  ArrayListConstIterator ret(*this);
  ++ (*this);
  return ret;
}

template<class T, unsigned int blocksize>
inline T const & ArrayListConstIterator<T, blocksize>::operator*() const  // dereferencing
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(!is_end());
  assert(index < nodeptr->used_size);
#endif
  return *(dataptr + index);
}

template<class T, unsigned int blocksize>
inline T const * ArrayListConstIterator<T, blocksize>::operator->() const
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(!is_end());
  assert(index < nodeptr->used_size);
#endif
  return (dataptr + index);
}


/**
  ArrayList
*/

template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::front() -> reference
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(next != nullptr);
  assert(next->used_size > 0);
  assert(num_blocks > 0);
  assert(total_size > 0);
  assert(static_cast<Block*>(next)->memory_buf != nullptr);
#endif
  return *(reinterpret_cast<T*>( static_cast<Block*>(next)->memory_buf) + (next->used_size) - 1);
}

template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::front() const -> const_reference
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(next != nullptr);
  assert(next->used_size > 0);
  assert(num_blocks > 0);
  assert(total_size > 0);
  assert(static_cast<Block*>(next)->memory_buf != nullptr);
#endif
  return *(reinterpret_cast<T const*>( (static_cast<Block*>(next)->memory_buf)) + (next->used_size) - 1);
}

template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::back() -> reference
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(prev != nullptr);
  assert(prev->used_size > 0);
  assert(num_blocks > 0);
  assert(total_size > 0);
  assert(static_cast<Block*>(prev)->memory_buf != nullptr);
#endif
  return *reinterpret_cast<T*>( static_cast<Block*>(prev)->memory_buf);
}

template<class T, unsigned int blocksize>
inline auto ArrayList<T,blocksize>::back() const -> const_reference
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(prev != nullptr);
  assert(prev->used_size > 0);
  assert(num_blocks > 0);
  assert(total_size > 0);
  assert(static_cast<Block*>(prev)->memory_buf != nullptr);
#endif
  return *reinterpret_cast<T const*>( static_cast<Block*>(prev)->memory_buf);
}

template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::cbegin() noexcept -> const_iterator
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(next != nullptr);
#endif
  return ArrayListConstIterator<T, blocksize>{next};
}

template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::cend() noexcept -> const_iterator
{
  return get_sentinel_iterator_const();
}

template<class T, unsigned int blocksize>
constexpr auto ArrayList<T, blocksize>::get_sentinel_iterator_const() noexcept -> const_iterator
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(used_size == 2*blocksize + 1);
#endif
  return ArrayListConstIterator<T, blocksize>{this, 0, nullptr};
}

template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::split_full_block_for_insert_before(const_iterator &pos) -> const_iterator
{
  Block & old_block = *static_cast<Block*>(pos.nodeptr);
  unsigned int const insertion_index = pos.index;
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(old_block.used_size == blocksize);
#endif
  // create new (empty) block preceding the current one.
  Block * const new_block = insert_block_between(old_block.prev, &old_block);
  static constexpr unsigned int new_size1 = (blocksize + 1)/2;  // new size of the old block
  static constexpr unsigned int new_size2 = blocksize + 1 - new_size1; // new size of the new block
  T * const old_mem = pos.dataptr;
  T * const new_mem = reinterpret_cast<T*>(new_block->memory_buf);
  if(insertion_index >= new_size1-1)
  {
    // old_block[0],...,old_block[new_size1 -1] can remain as they are.
    // old_block[newsize1], ... old_block[insertion_index] need to be moved into
    //   new_block[0] ,..., new_block[insertion_index - newsize1]
    // old_block[insertion_index+1],... need to be moved into
    //   new_block[insertion_index - newsize1 + 2]
    // Note that in the boundary case insertion_index == new_size1 - 1,
    // the insertion_index is the largest-index element in the old block, but the gap is in the new
    // block.
#ifdef USE_MEMMOVE
    // UB, but probably works...
    // TODO: At least add some static_asserts
    std::memcpy(new_mem, old_mem+new_size1, insertion_index-(new_size1-1) );
    std::memcpy(new_mem+insertion_index + 2 - new_size1, old_mem + insertion_index + 1, blocksize - insertion_index - 1 );
#else // ought to be equivalent
    for(unsigned int i=0; i < insertion_index - (new_size1 - 1); ++i )
    {
      // move by using placement new with move constructor and destroying the source.
      ::new(new_mem+i) T (std::move(old_mem[new_size1 + i]));
      old_mem[new_size1+i].~T();
    }
    for(unsigned int i=insertion_index + 1; i < blocksize; ++i)
    {
      ::new(new_mem + 1 - new_size1 + i) T (std::move(old_mem[i]));
      old_mem[i].~T();
    }
#endif
  }
  else // i.e insertion_index < new_size1 - 1
  {
    assert(false);
  }
  assert(false);
  return pos; // TODO: Change
}



} // end namespace GaussSieve


#endif // include guard
