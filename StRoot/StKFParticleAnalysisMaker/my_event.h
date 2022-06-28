#ifndef my_event_hh
#define my_event_hh
#include "TObject.h"
#include "KFParticle.h"

class my_event 
{
private:
    std::vector<KFParticle> particles;
public:
    my_event() {}
    my_event(std::vector<KFParticle> _particles);
    virtual ~my_event() {}
    std::vector<KFParticle> GetParticles();
    bool IsEmptyEvent() {return (particles.size() == 0);}
    void push_back(KFParticle _par);

ClassDef(my_event,1)
};

#endif
