#include "MixedBuffer.h"
#include "TRandom.h"

void MixedBuffer::Add(my_event _new_event, int cen, float vertexz)
{
    assert((cen >= 0 && cen <= 4) && "Invalid Centrality!");
    assert((fabs(vertexz <= 80)) && "Verte x out of range!");
    int vzindex = static_cast<int>(floor(vertexz/10.)+8.);
    if (vzindex == 16) vzindex--;
    if ((int)events[cen][vzindex].size() <  buffer_size) events[cen][vzindex].push_back(_new_event);
    else 
    {
        int rndm = gRandom->Integer(buffer_size);
        events[cen][vzindex][rndm] = _new_event;
    }

    assert((int)events[cen][vzindex].size() <= buffer_size && "Buffer overflow!");
}

void MixedBuffer::Add_FIFO(my_event _new_event, int cen, float vertexz)
{
    assert((cen >= 0 && cen <= 4) && "Invalid Centrality!");
    assert((fabs(vertexz <= 80)) && "Verte x out of range!");
    int vzindex = static_cast<int>(floor(vertexz/10.)+8.);
    if (vzindex == 16) vzindex--;
    if ((int)events[cen][vzindex].size() <  buffer_size) events[cen][vzindex].push_back(_new_event);
    else 
    {   
        for (int i = 0; i < buffer_size-1; ++i) events[cen][vzindex][i] = events[cen][vzindex][i+1];
        events[cen][vzindex][buffer_size-1] = _new_event;
    }

    assert((int)events[cen][vzindex].size() <= buffer_size && "Buffer overflow!");
}

void MixedBuffer::Add_Reservoir(my_event _new_event, int cen, float vertexz)
{
    assert((cen >= 0 && cen <= 4) && "Invalid Centrality!");
    assert((fabs(vertexz <= 80)) && "Verte x out of range!");
    int vzindex = static_cast<int>(floor(vertexz/10.)+8.);
    if (vzindex == 16) vzindex--;

    nEventsProcessed[cen][vzindex]++;
    if ((int)events[cen][vzindex].size() <  buffer_size) events[cen][vzindex].push_back(_new_event);
    else 
    {   
        if (gRandom->Rndm() < (buffer_size*1.0)/nEventsProcessed[cen][vzindex])
        {
            int index = gRandom->Integer(buffer_size);
            events[cen][vzindex][index] = _new_event;
        }
    }

    assert((int)events[cen][vzindex].size() <= buffer_size && "Buffer overflow!");
}

bool MixedBuffer::IsEmpty(int cen, float vertexz)
{   
    int vzindex = static_cast<int>(floor(vertexz/10.)+8.);
    if (vzindex == 16) vzindex--;
    return (events[cen][vzindex].size() == 0);
} 

int MixedBuffer::TotalStorage()
{
    int total = 0;
    for (int i = 0;i < 4; i++) for (int j = 0; j < 16; j++) total += events[i][j].size();
    return total;
}

my_event MixedBuffer::Sample(int cen, float vertexz)
{   
    int vzindex = static_cast<int>(floor(vertexz/10.)+8.);
    if (vzindex == 16) vzindex--;
    assert(events[cen][vzindex].size() != 0 && "Cannot sample from empty cell!");
    int rndm  = gRandom->Integer(events[cen][vzindex].size());
    return events[cen][vzindex][rndm];
}

std::vector<my_event> MixedBuffer::Sample_All(int cen, float vertexz)
{
    int vzindex = static_cast<int>(floor(vertexz/10.)+8.);
    if (vzindex == 16) vzindex--;
    assert(events[cen][vzindex].size() != 0 && "Cannot sample from empty cell!");
    return events[cen][vzindex];
}