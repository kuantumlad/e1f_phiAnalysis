#include "ParticleHistoFiller.h"
#include "TLorentzVector.h"
#include "TMath.h"

ParticleHistoFiller::ParticleHistoFiller(){


}

ParticleHistoFiller::ParticleHistoFiller(const char* temp_name){
  temp_hname = temp_name;
}

ParticleHistoFiller::~ParticleHistoFiller(){


}


void ParticleHistoFiller::initalizeAllHistos(){

  h_el_plotter.initializeHistos(Form("el_%s",temp_hname));
  h_pr_plotter.initializeHistos(Form("pr_%s",temp_hname));
  h_kp_plotter.initializeHistos(Form("kp_%s",temp_hname));
  h_km_plotter.initializeHistos(Form("km_%s",temp_hname));
  h_km0_plotter.initializeHistos(Form("km0_%s",temp_hname));
  h_km3_plotter.initializeHistos(Form("km3_%s",temp_hname));
  h_kaons0_plotter.InitKaonHistos(Form("km0_%s",temp_hname));
  h_kaons3_plotter.InitKaonHistos(Form("km3_%s",temp_hname));    

}

void ParticleHistoFiller::FillElectronHisto( TLorentzVector lv_temp ){

  h_el_plotter.h_p->Fill(lv_temp.P());
  h_el_plotter.h_theta->Fill(lv_temp.Theta() * (180.0/3.141592658535));
  h_el_plotter.h_phi->Fill(lv_temp.Phi()* (180.0/3.141592658535));

  h_el_plotter.h2_ptheta->Fill(lv_temp.P(), lv_temp.Theta() * (180.0/3.141592658535));
  h_el_plotter.h2_pphi->Fill(lv_temp.P(), lv_temp.Phi()* (180.0/3.141592658535));
  h_el_plotter.h2_thetaphi->Fill(lv_temp.Theta() * (180.0/3.141592658535), lv_temp.Phi()* (180.0/3.141592658535));
  
}

void ParticleHistoFiller::FillProtonHisto( TLorentzVector lv_temp ){

  h_pr_plotter.h_p->Fill(lv_temp.P());
  h_pr_plotter.h_theta->Fill(lv_temp.Theta() * (180.0/3.141592658535));
  h_pr_plotter.h_phi->Fill(lv_temp.Phi()* (180.0/3.141592658535));

  h_pr_plotter.h2_ptheta->Fill(lv_temp.P(), lv_temp.Theta() * (180.0/3.141592658535));
  h_pr_plotter.h2_pphi->Fill(lv_temp.P(), lv_temp.Phi()* (180.0/3.141592658535));
  h_pr_plotter.h2_thetaphi->Fill(lv_temp.Theta() * (180.0/3.141592658535), lv_temp.Phi()* (180.0/3.141592658535));
  
}

void ParticleHistoFiller::FillKaonPlusHisto( TLorentzVector lv_temp ){

  h_kp_plotter.h_p->Fill(lv_temp.P());
  h_kp_plotter.h_theta->Fill(lv_temp.Theta() * (180.0/3.141592658535));
  h_kp_plotter.h_phi->Fill(lv_temp.Phi()* (180.0/3.141592658535));

  h_kp_plotter.h2_ptheta->Fill(lv_temp.P(), lv_temp.Theta() * (180.0/3.141592658535));
  h_kp_plotter.h2_pphi->Fill(lv_temp.P(), lv_temp.Phi()* (180.0/3.141592658535));
  h_kp_plotter.h2_thetaphi->Fill(lv_temp.Theta() * (180.0/3.141592658535), lv_temp.Phi()* (180.0/3.141592658535));
  
}

void ParticleHistoFiller::FillKaonMinusHisto( TLorentzVector lv_temp ){

  h_km_plotter.h_p->Fill(lv_temp.P());
  h_km_plotter.h_theta->Fill(lv_temp.Theta() * (180.0/3.141592658535));
  h_km_plotter.h_phi->Fill(lv_temp.Phi()* (180.0/3.141592658535));

  h_km_plotter.h2_ptheta->Fill(lv_temp.P(), lv_temp.Theta() * (180.0/3.141592658535));
  h_km_plotter.h2_pphi->Fill(lv_temp.P(), lv_temp.Phi()* (180.0/3.141592658535));
  h_km_plotter.h2_thetaphi->Fill(lv_temp.Theta() * (180.0/3.141592658535), lv_temp.Phi()* (180.0/3.141592658535));  

}

void ParticleHistoFiller::FillKaonMinus0Histo( TLorentzVector lv_temp ){

  h_km0_plotter.h_p->Fill(lv_temp.P());
  h_km0_plotter.h_theta->Fill(lv_temp.Theta() * (180.0/3.141592658535));
  h_km0_plotter.h_phi->Fill(lv_temp.Phi()* (180.0/3.141592658535));

  h_km0_plotter.h2_ptheta->Fill(lv_temp.P(), lv_temp.Theta() * (180.0/3.141592658535));
  h_km0_plotter.h2_pphi->Fill(lv_temp.P(), lv_temp.Phi()* (180.0/3.141592658535));
  h_km0_plotter.h2_thetaphi->Fill(lv_temp.Theta() * (180.0/3.141592658535), lv_temp.Phi()* (180.0/3.141592658535));  

}

void ParticleHistoFiller::FillKaonMinus3Histo( TLorentzVector lv_temp ){

  h_km3_plotter.h_p->Fill(lv_temp.P());
  h_km3_plotter.h_theta->Fill(lv_temp.Theta() * (180.0/3.141592658535));
  h_km3_plotter.h_phi->Fill(lv_temp.Phi()* (180.0/3.141592658535));

  h_km3_plotter.h2_ptheta->Fill(lv_temp.P(), lv_temp.Theta() * (180.0/3.141592658535));
  h_km3_plotter.h2_pphi->Fill(lv_temp.P(), lv_temp.Phi()* (180.0/3.141592658535));
  h_km3_plotter.h2_thetaphi->Fill(lv_temp.Theta() * (180.0/3.141592658535), lv_temp.Phi()* (180.0/3.141592658535)); 

}


void ParticleHistoFiller::FillKaonPairHisto( std::vector<TLorentzVector> v_lv ){

  //MUST BE IN FOLLOWING ORDER
  TLorentzVector pr = v_lv[0];
  TLorentzVector kp = v_lv[1];
  TLorentzVector km = v_lv[2];
  
  double ang_prkm = pr.Vect().Angle(km.Vect());
  double ang_kpkm = kp.Vect().Angle(km.Vect());

  double pr_p = pr.P();
  double kp_p = kp.P();
  double km_p = km.P();
  

  //  h_kaons
  
  
}

void ParticleHistoFiller::FillKaonPair0Histo( std::vector<TLorentzVector> v_lv ){

  //MUST BE IN FOLLOWING ORDER
  TLorentzVector pr = v_lv[0];
  TLorentzVector kp = v_lv[1];
  TLorentzVector km = v_lv[2];
  
  double ang_prkm = pr.Vect().Angle(km.Vect()) * (180.0/TMath::Pi());
  double ang_kpkm = kp.Vect().Angle(km.Vect()) * (180.0/TMath::Pi());

  double pr_p = pr.P();
  double kp_p = kp.P();
  double km_p = km.P();
  double kp_km_p = kp_p;
  double prkm_p =  kp_p;
  //std::cout << " >> " << ang_kpkm << std::endl;
  h_kaons0_plotter.h_decayang->Fill( ang_kpkm );
  h_kaons0_plotter.h_decayangp->Fill( kp_km_p, ang_kpkm );
  h_kaons0_plotter.h_decayangp_prkm->Fill( prkm_p, ang_prkm );
  h_kaons0_plotter.h_decayang_prkm_kpkm->Fill( ang_kpkm, ang_prkm );
      
}

void ParticleHistoFiller::FillKaonPair3Histo( std::vector<TLorentzVector> v_lv ){

  //MUST BE IN FOLLOWING ORDER
  TLorentzVector pr = v_lv[0];
  TLorentzVector kp = v_lv[1];
  TLorentzVector km = v_lv[2];
  
  double ang_prkm = pr.Vect().Angle(km.Vect()) * (180.0/TMath::Pi());
  double ang_kpkm = kp.Vect().Angle(km.Vect()) * (180.0/TMath::Pi());

  double pr_p = pr.P();
  double kp_p = kp.P();
  double km_p = km.P();
  double kp_km_p = kp_p;
  double prkm_p = kp_p;
  
  h_kaons3_plotter.h_decayang->Fill( ang_kpkm );
  h_kaons3_plotter.h_decayangp->Fill( kp_km_p, ang_kpkm );
  h_kaons3_plotter.h_decayangp_prkm->Fill( prkm_p, ang_prkm );
  h_kaons3_plotter.h_decayang_prkm_kpkm->Fill( ang_kpkm, ang_prkm );
  
  
}

void ParticleHistoFiller::PrintAllHisto(){

  h_el_plotter.DisplayHistos(Form("el_%s",temp_hname));
  h_pr_plotter.DisplayHistos(Form("pr_%s",temp_hname));
  h_kp_plotter.DisplayHistos(Form("kp_%s",temp_hname));
  h_km_plotter.DisplayHistos(Form("km_%s",temp_hname));
  h_km0_plotter.DisplayHistos(Form("km0_%s",temp_hname));
  h_km3_plotter.DisplayHistos(Form("km3_%s",temp_hname));

  h_kaons0_plotter.DisplayHistos(Form("km0_%s",temp_hname));
  h_kaons3_plotter.DisplayHistos(Form("km3_%s",temp_hname));

}



