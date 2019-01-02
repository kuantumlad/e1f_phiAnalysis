#ifndef calcxparticle_hh
#define calcxparticle_hh

TLorentzVector CalcXParticle(std::vector< TLorentzVector > v_lv ){

  TLorentzVector lv_X = v_lv[0] + v_lv[1];
  for( int i = 2; i < v_lv.size(); i++ ){
    lv_X -= v_lv[i];
  }
  return lv_X;

};
#endif
