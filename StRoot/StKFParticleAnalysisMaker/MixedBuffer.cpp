#include "MixedBuffer.h"
#include "TRandom.h"

void MixedBuffer::Add(my_event _new_event, int cen, float vertexz)
{
    assert((cen >= 1 && cen <= 9) && "Invalid Centrality!");
    assert((fabs(vertexz <= 80)) && "Verte x out of range!");
    int vzindex = static_cast<int>(floor(vertexz/2.)+40.);
    if (vzindex == 80) vzindex--;
    if ((int)events[cen-1][vzindex].size() <  buffer_size) events[cen-1][vzindex].push_back(_new_event);
    else 
    {
        int rndm = gRandom->Integer(buffer_size);
        events[cen-1][vzindex][rndm] = _new_event;
    }

    assert((int)events[cen-1][vzindex].size() <= buffer_size && "Buffer overflow!");
}

void MixedBuffer::Add_FIFO(my_event _new_event, int cen, float vertexz)
{
    assert((cen >= 1 && cen <= 9) && "Invalid Centrality!");
    assert((fabs(vertexz <= 80)) && "Verte x out of range!");
    int vzindex = static_cast<int>(floor(vertexz/2.)+40.);
    if (vzindex == 80) vzindex--;
    if ((int)events[cen-1][vzindex].size() <  buffer_size) events[cen-1][vzindex].push_back(_new_event);
    else 
    {   
        for (int i = 0; i < buffer_size-1; ++i) events[cen-1][vzindex][i] = events[cen-1][vzindex][i+1];
        events[cen-1][vzindex][buffer_size-1] = _new_event;
    }

    assert((int)events[cen-1][vzindex].size() <= buffer_size && "Buffer overflow!");
}

void MixedBuffer::Add_Reservoir(my_event _new_event, int cen, float vertexz, int eventsprocessed)
{
    assert((cen >= 1 && cen <= 9) && "Invalid Centrality!");
    assert((fabs(vertexz <= 80)) && "Verte x out of range!");
    int vzindex = static_cast<int>(floor(vertexz/2.)+40.);
    if (vzindex == 80) vzindex--;
    if ((int)events[cen-1][vzindex].size() <  buffer_size) events[cen-1][vzindex].push_back(_new_event);
    else 
    {   
        if (gRandom->Rndm() < (buffer_size*1.0)/eventsprocessed)
        {
            int index = gRandom->Integer(buffer_size);
            events[cen-1][vzindex][index] = _new_event;
        }
    }

    assert((int)events[cen-1][vzindex].size() <= buffer_size && "Buffer overflow!");
}

bool MixedBuffer::IsEmpty(int cen, float vertexz)
{   
    int vzindex = static_cast<int>(floor(vertexz/2.)+40.);
    if (vzindex == 80) vzindex--;
    return (events[cen-1][vzindex].size() == 0);
} 

my_event MixedBuffer::Sample(int cen, float vertexz)
{   
    int vzindex = static_cast<int>(floor(vertexz/2.)+40.);
    if (vzindex == 80) vzindex--;
    assert(events[cen-1][vzindex].size() != 0 && "Cannot sample from empty cell!");
    int rndm  = gRandom->Integer(events[cen-1][vzindex].size());
    return events[cen-1][vzindex][rndm];
}

std::vector<my_event> MixedBuffer::Sample_All(int cen, float vertexz)
{
    int vzindex = static_cast<int>(floor(vertexz/2.)+40.);
    if (vzindex == 80) vzindex--;
    assert(events[cen-1][vzindex].size() != 0 && "Cannot sample from empty cell!");
    return events[cen-1][vzindex];
}