#include "my_event.h"
#include "TObject.h"

my_event::my_event(std::vector<KFParticle> _particles)
{
    particles = _particles;
}

std::vector<KFParticle> my_event::GetParticles() {return particles; }

void my_event::push_back(KFParticle _par) { particles.push_back(_par); }