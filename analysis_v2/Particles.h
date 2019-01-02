#ifndef particles_hh
#define particles_hh

#include <vector>

class Particles{

 public:
  Particles();
  virtual ~Particles();

  void SetEventParticles(std::vector<TLorentzVector>);
  
 protected:
  std::vector<TLorentzVector> event_particles;


  
  

};

#endif
