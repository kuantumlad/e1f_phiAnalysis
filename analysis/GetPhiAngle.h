#ifndef getphiangle_h
#define getphiangle_h


double GetPhiAngle( TLorentzVector lv_ebeam, TLorentzVector lv_target, TLorentzVector lv_el, TLorentzVector lv_pr, TLorentzVector lv_kp, TLorentzVector lv_km, TLorentzVector lv_phi ){
  
  TLorentzVector lv_vph;
  TLorentzVector lv_b_el(lv_el), lv_b_pr(lv_pr), lv_b_kp(lv_kp), lv_b_km(lv_km), lv_b_phi(lv_phi);;
  lv_vph = lv_ebeam - lv_el;
  TLorentzVector lv_wsys = -(lv_vph + lv_target);
  TVector3 v_boost = lv_wsys.BoostVector();
  lv_vph.Boost(v_boost);
  lv_b_el.Boost(v_boost);
  lv_b_pr.Boost(v_boost);
  lv_b_phi.Boost(v_boost);
  lv_kp.Boost(v_boost);
  lv_km.Boost(v_boost);
    
  TVector3 v_boost_vph = lv_vph.BoostVector();
  TVector3 v_boost_pr = lv_pr.BoostVector();
  TVector3 v_boost_phi = (lv_kp + lv_km).BoostVector();
    
  TVector3 v_lepton = (lv_vph.Vect()).Cross(lv_b_el.Vect());
  TVector3 v_hadron = (lv_vph.Vect()).Cross(lv_b_phi.Vect());
    
  TVector3 n_lepton = v_lepton.Unit();
  TVector3 n_hadron = v_hadron.Unit();
  double c0 = v_lepton.Dot(lv_b_phi.Vect());
  double n0 = n_lepton.Dot(n_hadron);
  double n1 = n_lepton.Mag();
  double n2 = n_hadron.Mag();
  double n_ang = (180.0/3.14159265)* (c0/TMath::Abs(c0))* TMath::ACos(n0/(n1*n2));

  return n_ang;
  
}

#endif
