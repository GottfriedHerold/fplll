#ifndef lattice_point_list_class_cpp
#define lattice_point_list_class_cpp

#include "PointList.h"
#include <assert.h>

template <class ZT>
void PointListMultiThreaded<ZT>::unlink(PointListIterator<ZT> const & pos, GarbageBin<ZT> &gb)
{
    {
        assert(!(pos.p->is_sentinel_node()));
        //std::lock_guard<std::mutex> writelock(mutex_currently_writing);
        mutex_currently_writing.lock();
        if(pos.p->is_plain_node()) //otherwise, already deleted.
            {
                pos.p->nodestatus=static_cast<int>(PointListMTNode<ZT>::StatusBit::is_to_be_deleted);
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





#endif
