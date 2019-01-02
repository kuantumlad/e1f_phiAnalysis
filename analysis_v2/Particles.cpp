#include "Particles.h"
#include "TLorentzVector.h"
#include <iostream>

Particles::Particles(){


}


Particles::~Particles(){

  

}


void Particles::SetEventParticles( std::vector<TLorentzVector> temp_event_particles){

  event_particles = temp_event_particles;

}

