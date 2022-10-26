#ifndef MixedBuffer_hh
#define MixedBuffer_hh
#include "./MixedBuffer.h"
#include "./my_event.h"
#include <cassert>

class MixedBuffer
{
private:
    std::vector<my_event> events[9][80]; 
    int nEventsProcessed[9][80];
    const int buffer_size;  
public:
    MixedBuffer():nEventsProcessed(), buffer_size(10) {}
    MixedBuffer(int _buffer_size):nEventsProcessed(), buffer_size(_buffer_size) {}
    virtual ~MixedBuffer() {}
    void Init() {for (int i = 0;i < 9; i++) for (int j = 0; j < 80; j++) {events[i][j].resize(0); nEventsProcessed[i][j]=0;}}
    void Add(my_event _new_event, int cen, float vertexz);
    void Add_FIFO(my_event _new_event, int cen, float vertexz);
    void Add_Reservoir(my_event _new_event, int cen, float vertexz);
    bool IsEmpty(int cen, float vertexz);
    my_event Sample(int cen, float vertexz);
    std::vector<my_event> Sample_All(int cen, float vertexz);

ClassDef(MixedBuffer,1)
};

#endif