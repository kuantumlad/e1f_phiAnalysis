#ifndef modparticle_hh
#define modparticle_hh

TLorentzVector ModParticle(TLorentzVector lv, double mod_mass ){

  TLorentzVector mod_lv(lv);
  mod_lv.SetPxPyPzE( lv.Px(), lv.Py(), lv.Pz(),  sqrt(lv.Vect().Mag2() + pow(mod_mass,2) ) );	       	

  return mod_lv;
};
#endif
