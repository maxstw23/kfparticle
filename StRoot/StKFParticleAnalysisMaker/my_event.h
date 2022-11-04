#ifndef my_event_hh
#define my_event_hh
#include "TObject.h"
#include "StPicoEvent/StPicoTrack.h"

class my_event 
{
private:
    std::vector<StPicoTrack*> particles;
public:
    my_event() {}
    my_event(std::vector<StPicoTrack*> _particles);
    virtual ~my_event() {}
    std::vector<StPicoTrack*> GetParticles();
    bool IsEmptyEvent() {return (particles.size() == 0);}
    void push_back(StPicoTrack* _par);

ClassDef(my_event,1)
};

#endif
