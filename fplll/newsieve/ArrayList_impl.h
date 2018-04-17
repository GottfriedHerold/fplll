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
    base_dataptr = reinterpret_cast<T*> (static_cast<Block*>(nodeptr)->memory_buf);
//    index   = nodeptr->used_size - 1;
    real_dataptr = base_dataptr + nodeptr->used_size - 1;
  }

  else // past-the end or before-begin iterator
  {
#ifdef DEBUG_SIEVE_ARRAYLIST
    assert(nodeptr->used_size == 2*blocksize + 1);
#endif
    base_dataptr = nullptr;
//    index   = 0;
    real_dataptr = nullptr;
  }
}

template<class T, unsigned int blocksize>
inline ArrayListConstIterator<T, blocksize>::ArrayListConstIterator(Meta* const new_nodeptr, unsigned int const new_index) noexcept
: // index(new_index),
  // dataptr(reinterpret_cast<T*>(static_cast< Block* >(new_nodeptr)->memory_buf)),
  // nodeptr(new_nodeptr)
  real_dataptr(reinterpret_cast<T*>(static_cast<Block*>(new_nodeptr)->memory_buf) + new_index),
  base_dataptr(reinterpret_cast<T*>(static_cast<Block*>(new_nodeptr)->memory_buf)),
  nodeptr(new_nodeptr)
{ }

template<class T, unsigned int blocksize>
constexpr inline ArrayListConstIterator<T, blocksize>::ArrayListConstIterator(Meta* const new_nodeptr, unsigned int const new_index, T * const new_base_dataptr) noexcept
: real_dataptr(new_base_dataptr + new_index), base_dataptr(new_base_dataptr), nodeptr(new_nodeptr)
{ }

template<class T, unsigned int blocksize>
inline ArrayListConstIterator<T,blocksize> &ArrayListConstIterator<T, blocksize>::operator++()
{
  if(real_dataptr != base_dataptr) // not at end of block.
  {
    --real_dataptr; // recall that index counts downwards
    return *this;
  }
  else // jump to next block
  {
    *this = ArrayListConstIterator{nodeptr->next};
//    if (nodeptr->next) { __builtin_prefetch(static_cast<Block*>(nodeptr->next)->memory_buf); }
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
  assert(get_index() < nodeptr->used_size);
#endif
  return *get_real_dataptr();
}

template<class T, unsigned int blocksize>
inline T const * ArrayListConstIterator<T, blocksize>::operator->() const
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(!is_end());
  assert(get_index() < get_nodeptr()->used_size);
#endif
  return get_real_dataptr();
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
  unsigned int const insertion_index = pos.get_index();
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(old_block->used_size == blocksize);
#endif
  // create new (empty) block preceding the current one.
  Block * const new_block = insert_block_between(old_block->prev, old_block);
  static constexpr unsigned int new_size1 = (blocksize + 1)/2;  // new size of the old block
  static constexpr unsigned int new_size2 = blocksize + 1 - new_size1; // new size of the new block
  T * const old_mem = pos.get_base_dataptr();
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
  unsigned int const insertion_index = pos.get_index();
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(old_block->used_size == blocksize);
#endif
  // create new (empty) block preceding the current one.
  Block * const new_block = insert_block_between(old_block->prev, old_block);
  static constexpr unsigned int new_size1 = (blocksize + 1)/2;  // new size of the old block
  static constexpr unsigned int new_size2 = blocksize + 1 - new_size1; // new size of the new block
  T * const old_mem = pos.get_base_dataptr();
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

/**
  Shift data around to create a gap in a given block.
*/

template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::create_gap_at(const_iterator &pos, unsigned int gap_index) -> const_iterator
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(pos.get_nodeptr() -> used_size < blocksize);
  assert(gap_index <= pos.nodeptr->used_size);
  assert(pos.get_base_dataptr() != nullptr);
#endif
#ifdef USE_MEMMOVE
  std::memmove(pos.get_base_dataptr() + gap_index + 1, pos.get_base_dataptr() + gap_index, sizeof(T) * (pos.get_nodeptr()->used_size - gap_index));
#else
  for(unsigned int i = pos.get_nodeptr()->used_size; i>=gap_index + 1; --i)
  {
    ::new(pos.get_base_dataptr() + i) T (pos.get_base_dataptr()[i-1]);
    pos.get_base_dataptr()[i-1].~T();
  }
#endif
  if(pos.get_index() >= gap_index)
  {
    // TODO: Encapsulation!
    ++pos.real_dataptr;
  }
  ++total_size;
  ++pos.nodeptr->used_size;
  return ArrayListConstIterator<T, blocksize>{pos.get_nodeptr(), gap_index, pos.get_base_dataptr()};
}

/**
  Removes an empty block
*/

template<class T, unsigned int blocksize>
inline void ArrayList<T,blocksize>::remove_empty_block(Block* block)
{
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(block->used_size == 0);
#endif
  block->prev->next = block->next;
  block->next->prev = block->prev;
  delete block;
  --num_blocks;
}

template<class T, unsigned int blocksize>
template<class... Args>
inline auto ArrayList<T,blocksize>::emplace_before(const_iterator &pos, Args&& ...args) -> const_iterator
{
  // Case 1 : There is still space inside the current node.
  // Note that this condition implies that we are not at the sentinel nodes, because
  // then we would have used_size > blocksize
  if(pos.get_nodeptr()->used_size < blocksize)
  {
//    std::cout << "Normal Insert";
#ifdef DEBUG_SIEVE_ARRAYLIST
    assert(pos.get_real_dataptr() != nullptr);
    assert(pos.get_base_dataptr() == reinterpret_cast<T*>(static_cast<Block*>(pos.get_nodeptr())->memory_buf) );
#endif
    auto retval = create_gap_at(pos, pos.get_index() + 1);
    ::new(retval.get_real_dataptr()) T (std::forward<Args>(args)...);
    return retval;
  }
  else if(pos.get_nodeptr()->used_size == blocksize) // full block
  {
//    std::cout << "Full block Insert";
    auto retval = split_full_block_for_insert_before(pos);
    ::new(retval.get_real_dataptr()) T (std::forward<Args>(args)...);
    return retval;
  }
  else // sentinel block, i.e. we insert at the END of the list
  {
    if(total_size == 0) // special case: yet-empty-list
    {
//      std::cout << "Initial Insert";
      // Checks that we really have an empty list.
      assert(pos.nodeptr == this);
      assert(next == this);
      assert(next == this);
      assert(num_blocks==0);

      Block * new_block = insert_block_between(this, this);
      ::new(new_block->memory_buf) T (std::forward<Args>(args)...);
      new_block->used_size = 1;
      total_size = 1;
      assert(num_blocks == 1); // set by insert_block_between
      // pos is still valid. We return the nemly inserted element
      return const_iterator{new_block, 0};
    }
    else // we move to the last block and insert there.
    {
//      std::cout << "Backed Insert";
      // const_iterator{pos.nodeptr->prev, 0} is an iterator to the last actual element.
      // Note that validity transformations to pos are unneeded.
      auto it = const_iterator{pos.nodeptr -> prev, 0}; // it may be modified by the call below
      return emplace_after(it, std::forward<Args>(args)... );
    }
  }
}

template<class T, unsigned int blocksize>
template<class... Args>
inline auto ArrayList<T,blocksize>::emplace_after(const_iterator &pos, Args&& ...args) -> const_iterator
{
  // Case 1: There is still free space in the current node. Just use that one.
  // (this implies that pos is no before-begin iterator)
  if(pos.get_nodeptr()->used_size < blocksize)
  {
//    std::cout << "BNormal Insert ";
#ifdef DEBUG_SIEVE_ARRAYLIST
    assert(pos.get_real_dataptr() != nullptr);
    assert(pos.get_real_dataptr() == reinterpret_cast<T*>(static_cast<Block*>(pos.nodeptr)->memory_buf) );
#endif
    auto retval = create_gap_at(pos, pos.get_index() );
    ::new(retval.get_real_dataptr()) T (std::forward<Args>(args)...);
    return retval;
  }
  else if(pos.get_nodeptr()->used_size == blocksize) // full block
  {
//    std::cout << "BFull Insert";
    auto retval = split_full_block_for_insert_after(pos);
    ::new(retval.get_real_dataptr() ) T (std::forward<Args>(args)...);
    return retval;
  }
  else // sentinel block, i.e. we insert at the END of the list
  {
    if(total_size == 0) // special case: yet-empty-list
    {
//      std::cout << "BInitialInsert";
      // Checks that we really have an empty list.
      assert(pos.get_nodeptr() == this);
      assert(next == this);
      assert(next == this);
      assert(num_blocks==0);

      Block * new_block = insert_block_between(this, this);
      ::new(new_block->memory_buf) T (std::forward<Args>(args)...);
      new_block->used_size = 1;
      total_size = 1;
      assert(num_blocks == 1); // set by insert_block_between
      // pos is still valid. We return the nemly inserted element
      return const_iterator{new_block, 0};
    }
    else // we move to the first block and insert there.
    {
//      std::cout << "ToBefore Insert";
      // const_iterator{pos.nodeptr->prev} is an iterator to the first actual element.
      // Note that validity transformations to pos are unneeded.
      auto it = const_iterator{pos.get_nodeptr() -> next};
      return emplace_before(it, std::forward<Args>(args)... );
    }
  }
}

template<class T, unsigned int blocksize>
inline auto ArrayList<T, blocksize>::erase(const_iterator &&pos) -> const_iterator
{
  Meta * const current_node = pos.get_nodeptr();
#ifdef DEBUG_SIEVE_ARRAYLIST
  assert(!(pos.is_end()));
  assert(pos.current_node->used_size > 0);
#endif
  --current_node->used_size;
#ifdef USE_MEMMOVE
  pos->~T();
  // +1 is performed on T* before converting to void*. Note that used_size was already decreased.
  std::memmove(pos.get_real_dataptr(),pos.get_real_dataptr() + 1,
               sizeof(T) * current_node->used_size - pos.get_index() );
#else
  // p is the target address for the move, which needs to range from the insertion position to
  // current_block[used_size - 1] (the new value of used_size)
  for(T *p = pos.get_real_dataptr(); p < pos.get_base_dataptr() + current_node -> used_size; ++p)
  {
    p->~T();
    ::new(p) T (std::move(*(p+1)));
  }
#endif
  /**
    WARNING: Fragile code!
  */
  // pos now refers to the element pseudo-_preceding_ the erasure position:
  // "pseudo", because if pos was the first element of the block (i.e. highest index),
  // pos is now invalid (refering to current_block[new_used_size] rather than the prev. block)
  // In particular, if the block is now empty, it is refering to current_block[0].
  // Still, ++pos will do the right thing.
  ++pos;
  // Note that current_node still points to the old node, in case ++pos switched block.
  if(current_node->used_size == 0)
  {
    remove_empty_block(static_cast<Block*>(current_node));
  }
  return pos;
}

} // end namespace GaussSieve


#endif // include guard
