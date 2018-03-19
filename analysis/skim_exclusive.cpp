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
#include <map>

#include "ExclusiveEvents.C"


void skim_exclusive( ){

  std::cout << ">> STARTING E1F PHI ANALYSIS " << std::endl;

  TChain *fchain = new TChain("events");
  fchain->Add("/Users/bclary/Documents/work/thesis_work/e1f_phiAnalysis/clas-phi/data/events.root");
  ExclusiveEvents *ex_event = new ExclusiveEvents(fchain);
  Long64_t nentries = fchain->GetEntries();
  
  TH1D *h_phy0_mm2 = new TH1D("h_phy0_mm2","h_phy0_mm2",200, 0.3, 0.7);
  TH1D *h_phy3_mm2 = new TH1D("h_phy3_mm2","h_phy3_mm2",200, 0.3, 0.7);
  TH1D *h_phyall_mm2 = new TH1D("h_phyall_mm2","h_phyall_mm2",200, 0.3, 0.7);
  TH1D *h_phy0_imkk = new TH1D("h_phy0_imkk","h_phy0_imkk", 200, 0.9, 1.1);
  TH1D *h_phy3_imkk = new TH1D("h_phy3_imkk","h_phy3_imkk", 200, 0.9, 1.1);
  TH1D *h_phy_ekX = new TH1D("h_phy_ekX","h_phy_ekX",200, 0.5, 2.0);
  TH1D *h_phy_t =  new TH1D("h_phy_t","h_phy_t",200, 0.0, 1.5);
  TH1D *h_phy_angphi = new TH1D("h_phy_angphi","h_phy_angphi", 100, -180.0, 180.0 );
  TH2D *h_phy03_imkk_kpr = new TH2D("h_phy03_imkk_kpr","h_phy03_imkk_kpr",300,0.7,2.0, 300, 1.4, 3.0); 
  
  TH1D *h_phy_ang_kk = new TH1D("h_phy_ang_kk","h_phy_ang_kk", 200, 0.0, 100.0 );
  TH1D *h_phy3_ang_kk = new TH1D("h_phy3_ang_kk","h_phy3_ang_kk", 200, 0.0, 100.0 );
  TH2D *h_phy_epkXIM = new TH2D("h_phy_epkXIM","h_phy_epkXIM",300, 0.8, 1.3, 300, 0.3, 0.7);
  
  TH1D *h_alpha_bin_pr = new TH1D("h_alpha_bin_pr","h_alpha_bin_pr",100, 0.0, 1.0);
  TH1D *h_alpha_bin_kp = new TH1D("h_alpha_bin_kp","h_alpha_bin_kp",100, 0.0, 1.0);
  TH1D *h_alpha_bin_km = new TH1D("h_alpha_bin_km","h_alpha_bin_km",100, 0.0, 1.0);

  TH1D *h_asy_phi_P = new TH1D("h_asy_phi_P","h_asy_phi_P",10, -180.0, 180.0 );
  TH1D *h_asy_phi_N = new TH1D("h_asy_phi_N","h_asy_phi_N",10, -180.0, 180.0 );
  
  TH1D *h_asy_delta_P = new TH1D("h_asy_delta_P","h_asy_delta_P",10, -180.0, 180.0 );
  TH1D *h_asy_delta_N = new TH1D("h_asy_delta_N","h_asy_delta_N",10, -180.0, 180.0 );  
  
  TLorentzVector lv_ebeam;
  lv_ebeam.SetPxPyPzE( 0, 0, 5.498, 5.498);
  TLorentzVector lv_target;
  lv_target.SetPxPyPzE( 0, 0, 0, 0.938);
  
  for( Long64_t i = 0; i < nentries; i++ ){
    fchain->GetEntry(i);

    int helicity = ((ex_event->helicity));
    int topology = ((ex_event->topology));
    
    TLorentzVector lv_el = (*(ex_event->electron));
    TLorentzVector lv_pr = (*(ex_event->proton));
    TLorentzVector lv_kp = (*(ex_event->kaon_pos));    
    TLorentzVector lv_km = (*(ex_event->kaon_neg));
    TLorentzVector lv_all_phi = lv_kp + lv_km;
    TLorentzVector lv_kpr = lv_pr + lv_kp;
    
    TLorentzVector lv_ekX = lv_ebeam + lv_target - lv_el - lv_kp;
    
    double t_pr = (lv_pr - lv_target).M2();

    TLorentzVector lv_vph;
    TLorentzVector lv_b_el(lv_el), lv_b_pr(lv_pr), lv_b_kp(lv_kp), lv_b_km(lv_km), lv_b_phi(lv_all_phi);;
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
    
    if( topology == 0 ){
      TLorentzVector lv_phi = lv_kp + lv_km;
      double ang_kk = (180.0/3.14159265) * (lv_kp.Vect().Angle(lv_km.Vect()) );
      h_phy_ang_kk->Fill(ang_kk);
      h_phy0_mm2->Fill(lv_km.M());
      h_phy_epkXIM->Fill( lv_phi.M(), lv_km.M() );
      h_phy_ekX->Fill( lv_ekX.M() );

      if( lv_ekX.M() >= 1.53 && lv_ekX.M() <= 1.57 ){
 	
	//if(
	if(helicity > 0){
	  h_asy_delta_P->Fill(n_ang);
	}
	else{
	  h_asy_delta_N->Fill(n_ang);
	}
      }
      
    }
    else if( topology == 3 ){
      TLorentzVector lv_phi = lv_kp + lv_km;
      double ang_kk = (180.0/3.14159265) * (lv_kp.Vect().Angle(lv_km.Vect()) );
      h_phy3_ang_kk->Fill( ang_kk );
      h_phy3_imkk->Fill(lv_phi.M());
      h_phy3_mm2->Fill(lv_km.M());
      h_phy_angphi->Fill(n_ang);

      // if( helicity > 0 ){
      //	h_asy_phi_P->Fill(n_ang);
      //}
      //else if ( helicity < 0 ) {
      //	h_asy_phi_N->Fill(n_ang);
      //}
      
    }
    TLorentzVector lv_phi_all = lv_kp + lv_km;
    h_phy03_imkk_kpr->Fill( lv_phi_all.M(), lv_ekX.M() );

    h_phy_t->Fill(-t_pr);
    
  }

  TH1D *h_asy_phi_top = new TH1D("h_asy_phi_top","h_asy_phi_top",10,-180.0, 180.0);// = h_asy_phi_P;
  TH1D *h_asy_phi_bot = new TH1D("h_asy_phi_bot","h_asy_phi_bot",10,-180.0, 180.0);;// = h_asy_phi_P;
  
  h_asy_phi_top->Add(h_asy_phi_P, h_asy_phi_N, -1.0);
  h_asy_phi_bot->Add(h_asy_phi_P, h_asy_phi_N, 1.0);
  
  TH1D *h_asy_phi_final =  new TH1D("h_asy_phi_final","h_asy_phi_final",10,-180.0, 180.0);
  h_asy_phi_final->Divide(h_asy_phi_top, h_asy_phi_bot, 1.0, 1.0);
  
  TGraph *g_asy_phi_final = new TGraph(h_asy_phi_final);
  g_asy_phi_final->SetName("g_asy_phi_final");
  g_asy_phi_final->SetTitle("CRUDE PHI BSA");


  TH1D *h_asy_delta_top = new TH1D("h_asy_delta_top","h_asy_phi_top",10,-180.0, 180.0);// = h_asy_phi_P;
  TH1D *h_asy_delta_bot = new TH1D("h_asy_delta_bot","h_asy_phi_bot",10,-180.0, 180.0);;// = h_asy_phi_P;
  
  h_asy_delta_top->Add(h_asy_delta_P, h_asy_delta_N, -1.0);
  h_asy_delta_bot->Add(h_asy_delta_P, h_asy_delta_N, 1.0);
  
  TH1D *h_asy_delta_final =  new TH1D("h_asy_delta_final","h_asy_delta_final",10,-180.0, 180.0);
  h_asy_delta_final->Divide(h_asy_delta_top, h_asy_delta_bot, 1.0, 1.0);
  
  TGraph *g_asy_delta_final = new TGraph(h_asy_delta_final);
  g_asy_delta_final->SetName("g_asy_delta_final");
  g_asy_delta_final->SetTitle("#Delta^{*} (1520) BSA");

  
  TCanvas *c_mm2 = new TCanvas("c_mm2","c_mm2",900,900);
  //c_mm2->Divide(1,1);
  //c_mm2->cd(1);
  h_phy_ekX->GetXaxis()->SetTitle("M.M. ek^{+}X");
  h_phy_ekX->GetXaxis()->CenterTitle();
  h_phy_ekX->Draw();
  //c_mm2->cd(2);
  //h_phy3_mm2->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",1600,800);
  c2->Divide(2,1);
  c2->cd(1);
  h_phy0_mm2->GetXaxis()->SetTitle("M.M. epk^{+}X [GeV]");
  h_phy0_mm2->GetXaxis()->CenterTitle();
  h_phy0_mm2->Draw();
  c2->cd(2);
  h_phy_epkXIM->GetXaxis()->SetTitle("IM_{kk} [GeV]");
  h_phy_epkXIM->GetXaxis()->CenterTitle();
  h_phy_epkXIM->GetYaxis()->SetTitle("M.M. epk^{+}X [GeV]");
  h_phy_epkXIM->GetYaxis()->CenterTitle();
  h_phy_epkXIM->Draw("colz");
  
  
  TCanvas *c4 = new TCanvas("c4","c4",900,900);
  gPad->SetLogz();
  h_phy03_imkk_kpr->GetXaxis()->SetTitle("I.M_{kk} [GeV]");
  h_phy03_imkk_kpr->GetXaxis()->CenterTitle();
  h_phy03_imkk_kpr->GetYaxis()->SetTitle("M.M. ek^{+}X [GeV]");
  h_phy03_imkk_kpr->GetYaxis()->CenterTitle();
  h_phy03_imkk_kpr->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  g_asy_delta_final->GetXaxis()->SetTitle("#phi [deg]");
  g_asy_delta_final->GetXaxis()->CenterTitle();
  g_asy_delta_final->GetYaxis()->SetTitle("BSA");
  g_asy_delta_final->GetYaxis()->CenterTitle();
  g_asy_delta_final->SetMarkerStyle(23);
  g_asy_delta_final->SetMarkerSize(1);
  g_asy_delta_final->Draw("AP");

}



//double getPhiAngle( TLorentzVector temp_el, TLorentzVector temp_pr, TLorentzVector
