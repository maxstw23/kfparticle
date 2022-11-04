#include "my_event.h"
#include "TObject.h"

my_event::my_event(std::vector<StPicoTrack*> _particles)
{
    particles = _particles;
}

std::vector<StPicoTrack*> my_event::GetParticles() {return particles; }

void my_event::push_back(StPicoTrack* _par) { particles.push_back(_par); }