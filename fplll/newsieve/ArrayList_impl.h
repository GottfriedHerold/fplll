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

/**
  split full block into two, insert and return gap BEFORE pos.
  pos is adjusted to account for memory movement.
*/


template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::split_full_block_for_insert_before(const_iterator &pos) -> const_iterator
{
  Block * const old_block = static_cast<Block*>(pos.nodeptr);
  unsigned int const insertion_index = pos.index;
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(old_block->used_size == blocksize);
#endif
  // create new (empty) block preceding the current one.
  Block * const new_block = insert_block_between(old_block->prev, old_block);
  static constexpr unsigned int new_size1 = (blocksize + 1)/2;  // new size of the old block
  static constexpr unsigned int new_size2 = blocksize + 1 - new_size1; // new size of the new block
  T * const old_mem = pos.dataptr;
  T * const new_mem = reinterpret_cast<T*>(new_block->memory_buf);
  old_block -> used_size = new_size1;
  new_block -> used_size = new_size2;
  ++total_size;
  if(insertion_index >= new_size1-1) // gap will be in the new block
  {
    // old_block[0],...,old_block[new_size1 -1] can remain as they are.
    // old_block[newsize1], ... old_block[insertion_index] need to be moved into
    // ->new_block[0] ,..., new_block[insertion_index - newsize1]
    // old_block[insertion_index+1],..., old_block[blocksize-1] need to be moved into
    // ->new_block[insertion_index - newsize1 + 2], ..., new_block[blocksize+1 - newsize1 - 1]
    // Note that in the boundary case insertion_index == new_size1 - 1,
    // the insertion_index is the largest-index element in the old block, but the gap is in the new
    // block.
#ifdef USE_MEMMOVE
    // UB, but probably works...
    // TODO: At least add some static_asserts
    std::memcpy(new_mem, old_mem+new_size1, sizeof(T)*(insertion_index-(new_size1-1)) );
    std::memcpy(new_mem+insertion_index + 2 - new_size1, old_mem + insertion_index + 1,
                sizeof(T)*(blocksize - insertion_index - 1 ));
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
    if(insertion_index > new_size1 -1) // old_block[insertion_index] was moved into new_block
    {
      // adjust pos to refer to the same element as before the splitting
      pos = ArrayListConstIterator<T, blocksize>{new_block,insertion_index-new_size1};
    }  // else insertion_index == new_size1 - 1 and pos remains unchanged
    // return position of gap
    return ArrayListConstIterator<T, blocksize>{new_block,insertion_index - new_size1 + 1, new_mem};

  }
  else // i.e insertion_index < new_size1 - 1
  {
    // old_block[0],..., old_block[insertion_index] can stay
    // old_block[insertion_index+1], ..., old_block[new_size-2] need to be moved to
    // ->old_block[insertion_index+2],..., old_block[new_size1-1]
    // old_block[new_size1-1],...,old_block[blocksize-1] moves to
    // ->new_block[0], ..., new_block[blocksize + 1 - new_size1 - 1] == new_block[new_size2-1]
#ifdef USE_MEMMOVE
    std::memcpy (new_mem, old_mem + new_size1 - 1, sizeof(T) * new_size2 );
    std::memmove(old_mem + insertion_index + 2, old_mem+insertion_index + 1, sizeof(T) * (new_size1 - insertion_index - 2));
#else
    for(unsigned int i=0; i < new_size2; ++i)
    {
      ::new(new_mem + i) T (std::move(old_mem[new_size1 - 1 + i] ) );
      old_mem[new_size1 - 1 + i].~T();
    }
    for(unsigned int i=new_size1 -1; i > insertion_index + 1; --i) // decreasing order!
    {
      ::new(old_mem + i) T (std::move(old_mem [i -1]));
      old_mem[i-1].~T();
    }
#endif
    // pos remains unchanged.
    return ArrayListConstIterator<T, blocksize>{old_block, insertion_index + 1, old_mem};
  }
}

/**
  split full block into two, insert and return gap AFTER pos.
  pos is adjusted to account for memory movement.
*/


template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::split_full_block_for_insert_after(const_iterator &pos) -> const_iterator
{
  Block * const old_block = static_cast<Block*>(pos.nodeptr);
  unsigned int const insertion_index = pos.index;
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(old_block->used_size == blocksize);
#endif
  // create new (empty) block preceding the current one.
  Block * const new_block = insert_block_between(old_block->prev, old_block);
  static constexpr unsigned int new_size1 = (blocksize + 1)/2;  // new size of the old block
  static constexpr unsigned int new_size2 = blocksize + 1 - new_size1; // new size of the new block
  T * const old_mem = pos.dataptr;
  T * const new_mem = reinterpret_cast<T*>(new_block->memory_buf);
  old_block -> used_size = new_size1;
  new_block -> used_size = new_size2;
  ++total_size;
  if(insertion_index >= new_size1) // position of new gap is in the new block
  {
    // old_block[0],...,old_block[new_size1 -1] can remain as they are.
    // old_block[newsize1], ... old_block[insertion_index-1] need to be moved into
    // ->new_block[0] ,..., new_block[insertion_index - newsize1 -1]
    // old_block[insertion_index],..., old_block[blocksize-1] need to be moved into
    // ->new_block[insertion_index - newsize1 + 1], ..., new_block[blocksize+1 - newsize1 - 1]
    // Gap is created at new_block[insertion_index - newsize1]
    // pos needs to be adjusted
#ifdef USE_MEMMOVE
    std::memcpy(new_mem, old_mem+new_size1, sizeof(T)*(insertion_index-new_size1));
    std::memcpy(new_mem+insertion_index + 1 - new_size1, old_mem + insertion_index,
                sizeof(T)*(blocksize - insertion_index ));
#else // ought to be equivalent
    for(unsigned int i=0; i < insertion_index - new_size1; ++i )
    {
      ::new(new_mem+i) T (std::move(old_mem[new_size1 + i]));
      old_mem[new_size1+i].~T();
    }
    for(unsigned int i=insertion_index; i < blocksize; ++i)
    {
      ::new(new_mem - new_size1 + i + 1) T (std::move(old_mem[i]));
      old_mem[i].~T();
    }
#endif
    pos = ArrayListConstIterator<T, blocksize>{new_block, insertion_index - new_size1 + 1, new_mem};
    // return position of gap
    return ArrayListConstIterator<T, blocksize>{new_block,insertion_index - new_size1, new_mem};
  }
  else // i.e insertion_index < new_size1
  {
    // old_block[0],..., old_block[insertion_index - 1] can stay
    // old_block[insertion_index], ..., old_block[new_size -2]  move to
    // ->old_block[insertion_index+1],..., old_block[new_size1-1]
    // old_block[new_size1-1],...,old_block[blocksize-1] moves to
    // ->new_block[0], ..., new_block[blocksize + 1 - new_size1 - 1] == new_block[new_size2-1]
    // NOTE: In the limiting case insertion_index  = new_size1 - 1, insertion_index is moved into
    // the new block by the 3rd case above!
#ifdef USE_MEMMOVE
    std::memcpy (new_mem, old_mem + new_size1 - 1, sizeof(T) * new_size2 );
    std::memmove(old_mem + insertion_index + 1, old_mem+insertion_index, sizeof(T) * (new_size1 - insertion_index - 1));
#else
    for(unsigned int i = 0; i < new_size2; ++i)
    {
      ::new(new_mem + i) T (std::move(old_mem[new_size1 - 1 + i]));
      old_mem[new_size1 -1 + i].~T();
    }
    for(unsigned int i = new_size1 -1; i > insertion_index; --i) // decreasing order!
    {
      ::new(old_mem + i) T (std::move(old_mem [i -1]));
      old_mem[i-1].~T();
    }
#endif
    if(insertion_index == new_size1 - 1)
    {
      pos = ArrayListConstIterator<T, blocksize>{new_block, 0, new_mem};
    }
    else
    {
      pos = ArrayListConstIterator<T, blocksize>{old_block, insertion_index + 1, old_mem};
    }
    return ArrayListConstIterator<T, blocksize>{old_block, insertion_index, old_mem};
  }
}

template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::create_gap_at(const_iterator &pos, unsigned int gap_index) -> const_iterator
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(pos.nodeptr -> used_size < blocksize);
  assert(gap_index <= pos.nodeptr->used_size);
  assert(pos.dataptr != nullptr);
#endif
#ifdef USE_MEMMOVE
  std::memmove(pos.dataptr + gap_index + 1, pos.dataptr + gap_index, sizeof(T) * (pos.nodeptr->used_size - gap_index));
#else
  for(unsigned int i = pos.nodeptr->used_size; i>=gap_index + 1; --i)
  {
    ::new(pos.dataptr + i) T (pos.dataptr[i-1]);
    pos.dataptr[i-1].~T();
  }
#endif
  if(pos.index >= gap_index)
  {
    ++pos.index;
  }
  return ArrayListConstIterator<T, blocksize>{pos.nodeptr, gap_index, pos.dataptr};
}

} // end namespace GaussSieve


#endif // include guard
