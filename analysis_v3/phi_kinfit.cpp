#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
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
#include "TPad.h"

#include "ExclusiveEvents.C"
#include "CLAS_KinFitter/include/KinFitter.h"
#include "CLAS_KinFitter/include/res-errs_clas.h"
//#include "ApplyKinFit.h"

void phi_kinfit(){

  std::cout << ">> STARTING E1F PHI ANALYSIS " << std::endl;

  TChain *fchain = new TChain("events");
  fchain->Add("/Users/bclary/Documents/work/thesis_work/e1f_phiAnalysis/clas-phi/data/events.root");
  ExclusiveEvents *ex_event = new ExclusiveEvents(fchain);
  //Long64_t nentries = fchain->GetEntries();
  Long64_t nentries = 3000;

  double ang_min = 170.0;
  int ang_bins = 110;
  TH2D *h2_ang_decay = new TH2D("h2_ang_decay","h2_ang_decay",ang_bins, ang_min, 180.0, ang_bins, ang_min, 180.0 );
  TH2D *h2_lang_decay = new TH2D("h2_lang_decay","h2_lang_decay",50, 0.0, 60.0, 50, 0.0, 60.0 );
  TH2D *h2_lang_decay_sig = new TH2D("h2_lang_decay_sig","h2_lang_decay_sig",50, 0.0, 60.0, 50, 0.0, 60.0 );

  TH2D *h2_ang_imkk = new TH2D("h2_ang_imkk","h2_ang_imkk", 100, 0.90, 1.7, 100, 150.0, 200.0);
  TH2D *h2_ang_imkmpr = new TH2D("h2_ang_imkmpr","h2_ang_imkmpr", 100, 1.3, 1.9, 100, 150.0, 200.0);

  TH2D *h2_ang_imkk2 = new TH2D("h2_ang_imkk2","h2_ang_imkk2", 100, 0.90, 1.7, 100, 0.0, 60.0);
  TH2D *h2_ang_imkmpr2 = new TH2D("h2_ang_imkmpr2","h2_ang_imkmpr2", 100, 1.3, 1.9, 100, 0.0, 60.0);

  TH1D *h_ang_kk = new TH1D("h_ang_kk","h_ang_kk",200, 0.0, 200.0);
  TH1D *h_ang_kmpr = new TH1D("h_ang_kmpr","h_ang_kmlam",200, 0.0, 200.0);
  TH1D *h_ang_diff = new TH1D("h_ang_diff","h_ang_diff",200, -50.0, 50.0 );
  TH1D *h_mass = new TH1D("h_mass","h_mass",75, 0.95, 1.7);
  TH1D *h_mass_sig = new TH1D("h_mass_sig","h_mass_sig",75, 0.95, 1.7);

  TH1D *h_s = new TH1D("h_s","h_s",80,1.0,6.0);
  
  TH2D *h_mass_s = new TH2D("h_mass_s","h_mass_s", 100, 1.0, 2.5, 100, 1.0, 7.0);
  TH2D *h_diff_s = new TH2D("h_pfidd_s","h_pdiff_s", 100, 1.0, 7.0, 100, -0.01, 0.01);
  TH2D *h_diff_ang = new TH2D("h_pfidd_ang","h_pdiff_ang", 100, 0.0, 60.0 , 100, -0.01, 0.01);
  TH2D *h_mass_phi_lam_sig = new TH2D("h_mass_phi_lam_sig","h_mass_phi_lam_sig", 50, 0.95, 1.5, 50, 1.4, 1.8);
  TH2D *h_mass_phi_lam_bg = new TH2D("h_mass_phi_lam_bg","h_mass_phi_lam_bg", 50, 0.95, 1.5, 50, 1.4, 1.8);
  TH2D *h_diff_pang = new TH2D("h_diff_pang","h_diff_pang",100, -0.01, 0.01, 100, -0.01, 0.01 );

  TH1D *h_mass_phi_lamX = new TH1D("*h_mass_phi_lamX","*h_mass_phi_labX",50,0.95,1.5);
  TH1D *h_mass_phi_lamY = new TH1D("*h_mass_phi_lamY","*h_mass_phi_lamY",50,1.4,1.8);
  
  TH1D *h_phy_imkk = new TH1D("h_phy_imkk", "h_phy_imkk", 75, 0.95, 1.3);
  TH1D *h_phy_imkk_kf = new TH1D("h_phy_imkk_kf", "h_phy_imkk_kf", 75, 0.95, 1.3);
  TH1D *h_phy_imkk_kf_sig = new TH1D("h_phy_imkk_kf_sig", "h_phy_imkk_kf_sig", 75, 0.95, 1.3);
  TH1D *h_kf_conf = new TH1D("h_kf_conf","h_kf_conf",100,0.0,1.05);

  TH3D *h3_phy_mass_s = new TH3D("h3_phy_mass_s","h3_phy_mass_s", 50, 0.95, 1.5, 50, 1.4, 1.8, 50, 1.0, 6.0); 
  
  TLorentzVector lv_ebeam;
  lv_ebeam.SetPxPyPzE( 0, 0, 5.498, 5.498);
  TLorentzVector lv_target;
  lv_target.SetPxPyPzE( 0, 0, 0, 0.938);
  
  KinFitter *kinFit = nullptr;
  TString fname   = "/Users/bclary/Documents/work/thesis_work/e1f_phiAnalysis/clas-phi/analysis_v3/fittedEvents.root";

  TFile *fileout  = new TFile( fname,"RECREATE");
  TTree *treeout = new TTree("kinfit","Kin. Fitted Events");

  double hel = 0;
  TBranch* b_kinfit = treeout->Branch("KinFitEvents", "KinFitter", &kinFit);
  TBranch* b_hel    = treeout->Branch("hel", &hel);  
  
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

    double acut = 0.00;
    if( alpha_proton >= acut && alpha_kaon_pos >= acut && topology == 3){
      if( lv_ekpX.M() >= 1.540 || lv_ekpX.M() <= 1.500 ){
	h_phy_imkk->Fill(lv_phi.M());
      }
      
      
      //std::cout << " >> " << alpha_proton << std::endl;
      TLorentzVector *P_e = new TLorentzVector(lv_el);
      TLorentzVector *P_pr = new TLorentzVector(lv_pr);
      TLorentzVector *P_kp = new TLorentzVector(lv_kp);
      TLorentzVector *P_km = new TLorentzVector(lv_km);
      TLorentzVector *P_phi = new TLorentzVector(lv_phi);
      std::vector<Double_t> sigmas_e, sigmas_p, sigmas_kp, sigmas_km;    
      sigmas_e    = getResErr_DC( *P_e );
      sigmas_p    = getResErr_DC( *P_pr );
      sigmas_kp = getResErr_DC( *P_kp );
      sigmas_km = getResErr_DC( *P_km );

      KFParticle KF_km_temp;
      
      if( topology == 0 ){
	sigmas_e[0]*=0.8;
	sigmas_e[1]*=0.8;
	sigmas_e[2]*=0.8;
	
	sigmas_p[0]*=0.8;        
	sigmas_p[1]*=0.8;        
	sigmas_p[2]*=0.8;        
	
	sigmas_kp[0]*=0.8;
	sigmas_kp[1]*=0.8;
	sigmas_kp[2]*=0.8;
	P_km->SetVectM( P_km->Vect(), 0.493);      
	KF_km_temp = KFParticle(*P_km, { },"km", 1, 0); //,"km"); //IF NOT MEASURED PUT EMPTY VECTOR IN WITH { } FOR THE SIGMAS, AND ADD 1, 0 TO CONSTRUCTOR
      }
      else if( topology == 3 ){
	sigmas_e[0]*=0.5;
	sigmas_e[1]*=0.5;
	sigmas_e[2]*=0.75;
	
	sigmas_p[0]*=0.5;        
	sigmas_p[1]*=0.5;        
	sigmas_p[2]*=0.75;        
	
	sigmas_kp[0]*=0.55;
	sigmas_kp[1]*=0.55;
	sigmas_kp[2]*=0.75;
	
	sigmas_km[0]*=0.55;
	sigmas_km[1]*=0.55;
	sigmas_km[2]*=0.75;
	KF_km_temp = KFParticle(*P_km, sigmas_km, "km"); //,"km"); //IF NOT MEASURED PUT EMPTY VECTOR IN WITH { } FOR THE SIGMAS, AND ADD 1, 0 TO CONSTRUCTOR
      }
      
      double phimass = P_phi->M();
      P_phi->SetVectM( P_phi->Vect(), 1.019);
      
      KFParticle  KF_e     =  KFParticle(*P_e, sigmas_e,"e");
      KFParticle  KF_p     =  KFParticle(*P_pr, sigmas_p,"p");
      KFParticle  KF_kp     =  KFParticle(*P_kp, sigmas_kp,"kp");
      KFParticle KF_km = KF_km_temp; //KFParticle(*P_km, {},"km",1,0);
                 
      //KFParticle  KF_phi    =  KFParticle(P_phi, {}, "phi", 1, 0);      
      //KFParticle  KF_km    =  KFParticle(*P_km, sigmas_km,"km");
    
      std::vector<KFParticle> KFParticles;
      KFParticles.push_back(KF_e);
      KFParticles.push_back(KF_p);
      KFParticles.push_back(KF_kp);
      KFParticles.push_back(KF_km);
	
      kinFit = new KinFitter( KFParticles ); // Initiali Kin. Fit with list of particles
      kinFit->SetM_targ(0.938);              // Set the target to be the proton (Default is helium, M_targ = 3.74 GeV/c^2 )
      kinFit->SetBeam(5.498);
      kinFit->DoFitting(100);                // Do the fit 100 times

      double cl_kinfit = kinFit->confLevel;
      h_kf_conf->Fill(cl_kinfit);

      TLorentzVector lvkf_e = kinFit->KFParticlesOut[0].P;
      TLorentzVector lvkf_pr = kinFit->KFParticlesOut[1].P;
      TLorentzVector lvkf_kp = kinFit->KFParticlesOut[2].P;
      TLorentzVector lvkf_km = kinFit->KFParticlesOut[3].P;
      TLorentzVector lvkf_phi = lvkf_kp + lvkf_km;
      TLorentzVector lvkf_ekpX = lv_ebeam + lv_target - lvkf_e - lvkf_kp;
      
      


      //if( cl_kinfit >= 0.2 ){
	h_phy_imkk_kf_sig->Fill((lvkf_kp + lvkf_km).M());
	h_mass_phi_lam_sig->Fill( lvkf_phi.M(), lvkf_ekpX.M() );
	h_mass_phi_lamX->Fill( lvkf_phi.M() );
	h_mass_phi_lamY->Fill( lvkf_ekpX.M() );
	h2_lang_decay_sig->Fill( lvkf_kp.Vect().Angle(lvkf_km.Vect()) * 180.0/3.141592658, lvkf_pr.Vect().Angle(lvkf_km.Vect()) * 180.0/3.141592658 );
	h3_phy_mass_s->Fill(lvkf_phi.M(), lvkf_ekpX.M(), frameS.Mag2() );
		
	//}
	//else{
	h_mass_phi_lam_bg->Fill( lvkf_phi.M(), lvkf_ekpX.M() );
	h2_lang_decay->Fill( lvkf_kp.Vect().Angle(lvkf_km.Vect()) * 180.0/3.141592658, lvkf_pr.Vect().Angle(lvkf_km.Vect()) * 180.0/3.141592658 );
	h_phy_imkk_kf->Fill((lvkf_kp + lvkf_km).M());
	//}
      treeout->Fill();      
    }
    //    delete kinFit;      
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",900,900);
  c1->SetTitle("IMKK BACKGROUND");
  h_phy_imkk_kf->Draw();
  c1->SaveAs("h_imkk.pdf");
  
  TCanvas *c2a = new TCanvas("c2","c2",900,900);
  //gStyle->SetOptStat(0);
  //TPad *center_pad = new TPad("center_pad", "center_pad",0.0,0.0,0.6,0.6);
  //center_pad->Draw();

  //TPad *top_pad = new TPad("top_pad", "top_pad",0.0,0.55,0.6,1.0);
  //top_pad->Draw();  

  //TPad *right_pad = new TPad("right_pad", "right_pad",0.55,0.0,1.0,0.6);
  //right_pad->Draw();  
  
  //center_pad->cd();
  //gStyle->SetPalette(1);
  h_mass_phi_lam_sig->Draw("colz");

  //top_pad->cd();
  //TH1D *h_tempXa = new TH1D("h_tempXa","h_tempXa",50,0.95,1.5);
  //h_mass_phi_lam_sig->ProjectionX("h_tempXa",1,50);
  //h_tempXa->SetFillColorAlpha(kBlue+1,0.7);
  //h_tempXa->Draw("hist");
  //c2a->Update();

  //right_pad->cd();
  //h_mass_phi_lamY->SetFillColorAlpha(kBlue+2,0.7);
  //h_mass_phi_lamY->Draw();
  //c2a->Update();   
  c2a->SaveAs("h2_mass_phi_lab_sig.pdf");

  TCanvas *c2b = new TCanvas("c2b","c2b",900,900);
  //TH1D *h_temp2X = h_mass_phi_lam_sig->ProjectionX();
  //h_tempXa->Draw("hist");    
  h_mass_phi_lam_bg->Draw("colz");
  c2b->SaveAs("h2_mass_phi_lam_bg.pdf");

  
  
  
  TCanvas *c3 = new TCanvas("c3","c3",900,900);
  gPad->SetLogy();
  h_kf_conf->Draw();
  c3->SaveAs("h_kf_conf.pdf");

  TCanvas *c4 = new TCanvas("c4","c4",900,450);
  c4->Divide(2,1);
  c4->cd(1);
  h2_lang_decay->Draw("colz");
  c4->cd(2);
  h2_lang_decay_sig->Draw("colz");
  c4->SaveAs("h2_lab_ang_sb.pdf");

  TCanvas *c5 = new TCanvas("c5","c5",900,900);
  c5->Divide(1,1);
  c5->cd(1);
  h_phy_imkk_kf_sig->SetLineColor(kRed);
  h_phy_imkk_kf->Draw();
  h_phy_imkk_kf_sig->Draw("same");
  c5->SaveAs("h_kf_imkk_sb.pdf");

  TCanvas *c6 = new TCanvas("c6","c6",900,900);
  h_phy_imkk_kf->SetLineColor(kRed);
  h_phy_imkk->Draw();
  h_phy_imkk_kf->Draw("same");
  c6->SaveAs("h_imkk_comp.pdf");

  TCanvas *c7 = new TCanvas("c7","c7",900,900);
  h3_phy_mass_s->Draw("COL");
  
  fileout->Write();
  fileout->Close();
  delete fileout;
  std::cout << " DONE " << std::endl;
  
  return;
}
