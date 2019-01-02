#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TFile.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TString.h>
#include "TROOT.h"
#include "TBrowser.h"
#include <iostream>
#include <string>
#include "TChain.h"
#include "TRegexp.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TF1.h"
#include <map>
#include "ExclusiveEvents.C"


void ang_test(){

  std::cout << ">> STARTING E1F PHI ANALYSIS " << std::endl;

  TChain *fchain = new TChain("events");
  fchain->Add("/Users/bclary/Documents/work/thesis_work/e1f_phiAnalysis/clas-phi/data/events.root");
  ExclusiveEvents *ex_event = new ExclusiveEvents(fchain);
  Long64_t nentries = fchain->GetEntries();

  double ang_min = 170.0;
  int ang_bins = 110;
  TH2D *h2_ang_decay = new TH2D("h2_ang_decay","h2_ang_decay",ang_bins, ang_min, 180.0, ang_bins, ang_min, 180.0 );
  TH2D *h2_lang_decay = new TH2D("h2_lang_decay","h2_lang_decay",50, 0.0, 60.0, 50, 0.0, 60.0 );

  TH2D *h2_ang_imkk = new TH2D("h2_ang_imkk","h2_ang_imkk", 100, 0.90, 1.7, 100, 150.0, 200.0);
  TH2D *h2_ang_imkmpr = new TH2D("h2_ang_imkmpr","h2_ang_imkmpr", 100, 1.3, 1.9, 100, 150.0, 200.0);

  TH2D *h2_ang_imkk2 = new TH2D("h2_ang_imkk2","h2_ang_imkk2", 100, 0.90, 1.7, 100, 0.0, 60.0);
  TH2D *h2_ang_imkmpr2 = new TH2D("h2_ang_imkmpr2","h2_ang_imkmpr2", 100, 1.3, 1.9, 100, 0.0, 60.0);

  TH1D *h_ang_kk = new TH1D("h_ang_kk","h_ang_kk",200, 0.0, 200.0);
  TH1D *h_ang_kmpr = new TH1D("h_ang_kmpr","h_ang_kmlam",200, 0.0, 200.0);
  TH1D *h_ang_diff = new TH1D("h_ang_diff","h_ang_diff",200, -50.0, 50.0 );
  TH1D *h_mass = new TH1D("h_mass","h_mass",75, 0.95, 1.7);

  TH1D *h_s = new TH1D("h_s","h_s",80,1.0,6.0);
  
  TH2D *h_mass_s = new TH2D("h_mass_s","h_mass_s", 100, 1.0, 2.5, 100, 1.0, 7.0);
  TH2D *h_diff_s = new TH2D("h_pfidd_s","h_pdiff_s", 100, 1.0, 7.0, 100, -0.01, 0.01);
  TH2D *h_diff_ang = new TH2D("h_pfidd_ang","h_pdiff_ang", 100, 0.0, 60.0 , 100, -0.01, 0.01);
  TH2D *h_mass_phi_lam = new TH2D("h_mass_phi_lam","h_mass_phi_lam",100,0.9,1.5, 100, 1.1, 1.7);
  TH2D *h_diff_pang = new TH2D("h_diff_pang","h_diff_pang",100, -0.01, 0.01, 100, -0.01, 0.01 );
  
  TLorentzVector lv_ebeam;
  lv_ebeam.SetPxPyPzE( 0, 0, 5.498, 5.498);
  TLorentzVector lv_target;
  lv_target.SetPxPyPzE( 0, 0, 0, 0.938);

  std::cout << nentries << std::endl;
  
  
  for( Long64_t i = 0; i < nentries; i++ ){
    fchain->GetEntry(i);
    
    int helicity = ((ex_event->helicity));
    int topology = ((ex_event->topology));
    double alpha_proton = ((ex_event->alpha_proton));
    double alpha_kaon_pos = ((ex_event->alpha_kaon_pos));
    double alpha_kaon_neg = ((ex_event->alpha_kaon_neg));
    
    TLorentzVector lv_el = (*(ex_event->electron));       
    TLorentzVector lv_pr = (*(ex_event->proton));
    TLorentzVector lv_kp = (*(ex_event->kaon_pos));    
    TLorentzVector lv_km = (*(ex_event->kaon_neg));
    TLorentzVector lv_phi = lv_kp + lv_km;
    TLorentzVector lv_lambda = lv_pr + lv_km;
    TLorentzVector lv_epX = lv_ebeam + lv_target - lv_el - lv_pr;    
    TLorentzVector lv_ekpX = lv_ebeam + lv_target - lv_el - lv_kp;    
    
    TLorentzVector lv_virtualphoton = lv_ebeam - lv_el;
    TLorentzVector lv_b_el = lv_el;
    TLorentzVector lv_b_pr = lv_pr;
    TLorentzVector lv_b_kp = lv_kp;
    TLorentzVector lv_b_km = lv_km;
    TLorentzVector lv_b_km2 = lv_km;
    TLorentzVector lv_b_lam = lv_km + lv_pr;
    TLorentzVector lv_b_phi = lv_kp + lv_km;
    
    TLorentzVector frameW= -(lv_virtualphoton + lv_target);
    TLorentzVector frameS = (lv_el + lv_target);
    TLorentzVector framePHI = -(lv_ebeam + lv_target - lv_el - lv_pr);
    TLorentzVector frameLAMBDA = -(lv_ebeam + lv_target - lv_el - lv_kp);
      
    TVector3 boostW = -frameW.BoostVector();
    TVector3 boostS = frameS.BoostVector();
    TVector3 boostPHI = -framePHI.BoostVector();
    TVector3 boostLAMBDA = -frameLAMBDA.BoostVector();

    //std::cout << " >> " << boostPHI.Mag() << "  " << boostLAMBDA.Mag() << std::endl;
    lv_b_kp.Boost(boostPHI);
    lv_b_km.Boost(boostPHI);
        
    lv_b_pr.Boost(boostLAMBDA);
    lv_b_km2.Boost(boostLAMBDA);
    lv_b_lam.Boost(boostLAMBDA);
      
    double ang_kk = (lv_b_km.Vect()).Angle(lv_b_kp.Vect()) * 180.0/3.1415926538;
    double ang_kmpr = (lv_b_km2.Vect()).Angle(lv_b_pr.Vect()) * 180.0/3.1415926538;

    double lab_ang_kk = (lv_km.Vect()).Angle(lv_kp.Vect()) * 180.0/3.1415926538;
    double lab_ang_kmpr = (lv_km.Vect()).Angle(lv_pr.Vect()) * 180.0/3.1415926538;

    double kaon_diff = (lv_kp.Mag() - lv_km.Mag())/(lv_kp.Mag() + lv_km.Mag()) * 10;
    double ang_diff = (ang_kk - ang_kmpr )/(ang_kk + ang_kmpr) * 10;
    
    //    bool ang_cut_top = lab_ang_kmpr < 2.0 + 0.8*lab_ang_kk;
    //bool ang_cut_bot = lab_ang_kmpr > -7.0 + 0.8*lab_ang_kk;

    double acut = 0.4;
    //if( topology == 0 )
    //if ( alpha_kaon_pos > acut && alpha_proton > acut){
      
      h_ang_kk->Fill(ang_kk);
      h_ang_kmpr->Fill(ang_kmpr);
	
      //if( ang_cut_top && ang_cut_bot ) {//  && ang_kk <= 17.0 ){
      //if( ang_kk  > 175.0  || ang_kmpr < 175.0 ){

      //      if( ang_kk < 174.925  && ang_kmpr < 179.4){

	double s = (lv_el + lv_target).Mag2();
	h_s->Fill(s);      
	int sbin = h_s->GetBin(s);
	//std::cout << " << " << sbin << std::endl;
	if( sbin >= 0 ){
	  h_mass->Fill(lv_epX.M());
	  h2_ang_imkmpr->Fill(lv_lambda.M(), ang_kmpr);
	  h_mass_phi_lam->Fill( lv_phi.M(), lv_ekpX.M() ); 
	  h2_ang_decay->Fill(ang_kk, ang_kmpr);
	  
	  h2_ang_imkk->Fill(lv_phi.M(), ang_kk );
	  h2_lang_decay->Fill(lab_ang_kk, lab_ang_kmpr);	
	  h_mass_s->Fill( lv_ekpX.M(), s );
	  
	  h_diff_s->Fill(s, kaon_diff);
	  h_diff_ang->Fill(s, lab_ang_kk);
	  h_diff_pang->Fill(kaon_diff, ang_diff);	  
	}
	//      }
	//}

      //}
      //}
      //if( ang_cut_bot && ang_kmpr > 17.0 ){
      
      //if( ang_kmpr > 170 ){
	//}
      //}
	// }
  }

  TCanvas *c1 = new TCanvas("c1","c1",900,450);
  c1->Divide(2,1);
  c1->cd(1);
  h_ang_kk->Draw();
  c1->cd(2);
  //h_ang_kmpr->Add(h_ang_kk);
  h_ang_kmpr->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",900,450);
  gPad->SetLogz();
  c2->Divide(2,1);
  c2->cd(1);      
  TF1 *ang_cut = new TF1("ang_cut","3.0 + 0.8*x",0.0,200.0);
  h2_ang_decay->Draw("colz");
  ang_cut->Draw("same");
  c2->cd(2);
  h2_lang_decay->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","c3",900,450);
  c3->Divide(2,1);
  c3->cd(1);
  h2_ang_imkk->Draw("colz");
  c3->cd(2);
  h2_ang_imkmpr->Draw("colz");

  TCanvas *c4 = new TCanvas("c4","c4",900,450);
  c4->Divide(2,1);
  c4->cd(1);
  h_mass->Draw();
  c4->cd(2);
  h_mass_s->Draw("colz");

  TCanvas *c5 = new TCanvas("c5","c5",900,450);
  c5->Divide(2,1);
  c5->cd(1);
  h_s->Draw();
  c5->cd(2);
  h_mass_phi_lam->Draw("colz");

  TCanvas *c6 = new TCanvas("c6","c6",900,900);
  h_diff_s->Draw("colz");

  TCanvas *c7 = new TCanvas("c7","c7",900,450);
  c7->Divide(2,1);
  c7->cd(1);
  h_diff_ang->Draw("colz");
  c7->cd(2);
  h_diff_pang->Draw("colz");
   
  
  
}
