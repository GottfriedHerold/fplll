//continuation from PointList.h
//no include guards or anything needed.

template<class ET>
GaussList<ET,true,-1>::GaussList() :
    mutex_currently_writing(),
    start_sentinel_node (new Node),
    end_sentinel_node   (new Node)
{
        start_sentinel_node->next_node=end_sentinel_node;
        end_sentinel_node->prev_node=start_sentinel_node;
        start_sentinel_node->nodestatus=static_cast<int>(Node::StatusBit::is_first_node);
        end_sentinel_node->nodestatus  =static_cast<int>(Node::StatusBit::is_last_node);
}



template<class ET>
GaussList<ET,true,-1>::~GaussList() //called when only one thread is running
{
    //atomic_thread_fence(std::memory_order_acquire);   begin() already does that.
    Iterator next=cbegin();
    Iterator cur (start_sentinel_node);
    while(!next.is_end())
    {
        delete cur.p;
        cur=next;
        ++next;
    }
    delete cur.p;
    delete next.p;
}





















































//old implementation:


template<class DT>
MTListIterator<DT> MTListIterator<DT>::operator++(int) //postfix version
{
    auto tmp=p;
    ++(*this);
    return tmp;
};

template<class DT>
MTListIterator<DT>& MTListIterator<DT>::operator++()
{
    p=p->next_node.load(memory_order_acquire);
    assert(p!=nullptr);
    return *this;
}; //prefix version

/*
template<class DT>
typename ListMultiThreaded<DT>::Iterator ListMultiThreaded<DT>::cbefore_begin() const
{
    return start_sentinel_node; //no atomic load, because
}*/

/*
template<class DT>
typename ListMultiThreaded<DT>::Iterator ListMultiThreaded<DT>::cbegin() const
{
    return start_sentinel_node->next_node.load(std::memory_order_acquire);
};
*/

/*
template<class DT>
typename ListMultiThreaded<DT>::Iterator ListMultiThreaded<DT>::cend() const
{
    return end_sentinel_node;
};
*/

/*
template<class DT>
typename ListMultiThreaded<DT>::Iterator ListMultiThreaded<DT>::insert_before(Iterator const &pos, DT const &val)
{
    DT * const tmp = new DT(val);
    return enlist(pos,tmp);
};
*/

/*
template<class DT>
ListMultiThreaded<DT>::~ListMultiThreaded() //called when only one thread is running
{
    //atomic_thread_fence(std::memory_order_acquire);   begin() already does that.
    Iterator next=cbegin();
    Iterator cur (start_sentinel_node);
    while(!next.is_end())
    {
        delete cur.p;
        cur=next;
        ++next;
    }
    delete cur.p;
    delete next.p;
}
*/

/*
template<class DT>
ListMultiThreaded<DT>::ListMultiThreaded() : //called when only one thread is running
        mutex_currently_writing(),
        start_sentinel_node (new Node),
        end_sentinel_node   (new Node)
    {
        start_sentinel_node->next_node=end_sentinel_node;
        end_sentinel_node->prev_node=start_sentinel_node;
        start_sentinel_node->nodestatus=static_cast<int>(Node::StatusBit::is_first_node);
        end_sentinel_node->nodestatus  =static_cast<int>(Node::StatusBit::is_last_node);
    };
*/

/*
template <class DT>
typename ListMultiThreaded<DT>::Iterator ListMultiThreaded<DT>::enlist_before(MTListIterator<DT> const &pos, DT * const valref)
{
    Node* newnode = new Node;
    newnode->latpoint = valref;
    Node* nextgood =pos.p; //we work directly with the underlying pointer, not using the iterator, since we do not need atomic loads to traverse here.
    assert(nextgood!=nullptr);
    mutex_currently_writing.lock();
    while(nextgood->is_marked_for_deletion())
    {
        nextgood=nextgood->next_node;
    } //.load(memory_order_relaxed);}
    Node* preced=nextgood->prev_node;
    newnode->next_node=nextgood;
    newnode->prev_node=preced;
    nextgood->prev_node=newnode;
//until here, no other thread can observe our writes, since we did not publish the pointer.
    preced->next_node.store(newnode,memory_order_release);
//this one can be observed (even prior to releasing the lock), and we have to make sure that other non-mutex-protected threads that see this write also see the (non-atomic) writes to newnode.
    mutex_currently_writing.unlock();
    return static_cast<Iterator> (newnode);
}
*/

/*
template <class DT>
void ListMultiThreaded<DT>::unlink(Iterator const & pos, GarbageBin<DT> &gb)
{
    assert(pos.p!=nullptr);
    {//locked part
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
    }//end of locked part
    return;
}
*/

template<class DT>
void GarbageBin<DT>::empty_trash()
{
    while(!this->empty())
    {
        delete this->front();
        this->pop();
    }
};

