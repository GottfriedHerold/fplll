#ifndef lattice_point_list_class_cpp
#define lattice_point_list_class_cpp

#include "PointList.h"
#include <assert.h>

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
void ListMultiThreaded<DT>::enlist(Iterator const &pos, DT const * const &valref)
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
newnode.next_node=next_good.p;
newnode.prev_node=next_good.p->prev_node;
newnode.next_node->prev_node=newnode;
newnode.prev_node->next_node=newnode;
//does not work -- need atomics to prevent compiler from reorderings.
mutex_currently_writing.unlock();

}
#endif
