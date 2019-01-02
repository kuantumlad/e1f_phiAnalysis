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
//#include "CLAS_KinFitter/include/KinFitter.h"
//#include "CLAS_KinFitter/include/res-errs_clas.h"

void skim_exclusive( ){

  std::cout << ">> STARTING E1F PHI ANALYSIS " << std::endl;

  TChain *fchain = new TChain("events");
  fchain->Add("/Users/bclary/Documents/work/thesis_work/e1f_phiAnalysis/clas-phi/data/events.root");
  ExclusiveEvents *ex_event = new ExclusiveEvents(fchain);
  Long64_t nentries = fchain->GetEntries();
  
  TH1D *h_phy0_mm2 = new TH1D("h_phy0_mm2","h_phy0_mm2",200, 0.3, 0.7);//0.47, 0.53);//0.3, 0.7);
  TH1D *h_phy3_mm2 = new TH1D("h_phy3_mm2","h_phy3_mm2",200, 0.3, 0.7);
  TH1D *h_phyall_mm2 = new TH1D("h_phyall_mm2","h_phyall_mm2",200, 0.3, 0.7);
  TH1D *h_phy0_imkk = new TH1D("h_phy0_imkk","h_phy0_imkk", 200, 0.9, 1.1);
  TH1D *h_phy3_imkk = new TH1D("h_phy3_imkk","h_phy3_imkk", 200, 0.9, 1.1);
  TH1D *h_phy_ekX = new TH1D("h_phy_ekX","h_phy_ekX",200, 0.5, 2.0);
  TH1D *h_phy_ekX_ac = new TH1D("h_phy_ekX_ac","h_phy_ekX_ac",200, 0.5, 2.0);
  TH1D *h_phy_t =  new TH1D("h_phy_t","h_phy_t",200, 0.0, 1.5);
  TH1D *h_phy_angphi = new TH1D("h_phy_angphi","h_phy_angphi", 100, -180.0, 180.0 );
  TH2D *h_phy03_imkk_kpr = new TH2D("h_phy03_imkk_kpr","h_phy03_imkk_kpr",300,0.7,2.0, 300, 1.4, 3.0); 

  TH1D *h_phy03_imkk = new TH1D("h_phy03_imkk","h_phy03_imkk", 50, 0.95, 1.15);
  TH1D *h_phy03_imkk_kf = new TH1D("h_phy03_imkk_kf","h_phy03_imkk_kf",200, 0.95, 1.35);
  
  TH1D *h_phy_ang_kk = new TH1D("h_phy_ang_kk","h_phy_ang_kk", 200, 0.0, 100.0 );

  TH2D *h_phy_ang_decay = new TH2D("h_phy_ang_decay","h_phy_ang_decay", 500, -360.0, 360.0, 500, -360.5, 360.50 );

  TH1D *h_phy3_ang_kk = new TH1D("h_phy3_ang_kk","h_phy3_ang_kk", 200, 0.0, 100.0 );
  TH2D *h_phy_epkXIM = new TH2D("h_phy_epkXIM","h_phy_epkXIM",300, 0.8, 1.3, 300, 0.3, 0.7);
  
  TH1D *h_alpha_bin_pr = new TH1D("h_alpha_bin_pr","h_alpha_bin_pr",100, 0.0, 1.0);
  TH1D *h_alpha_bin_kp = new TH1D("h_alpha_bin_kp","h_alpha_bin_kp",100, 0.0, 1.0);
  TH1D *h_alpha_bin_km = new TH1D("h_alpha_bin_km","h_alpha_bin_km",100, 0.0, 1.0);

  TH1D *h_asy_phi_P = new TH1D("h_asy_phi_P","h_asy_phi_P",10, -180.0, 180.0 );
  TH1D *h_asy_phi_N = new TH1D("h_asy_phi_N","h_asy_phi_N",10, -180.0, 180.0 );
  
  TH1D *h_asy_delta_P = new TH1D("h_asy_delta_P","h_asy_delta_P",10, -180.0, 180.0 );
  TH1D *h_asy_delta_N = new TH1D("h_asy_delta_N","h_asy_delta_N",10, -180.0, 180.0 );

  TH1D *h_phy_kpi = new TH1D("h_phy_kpi","h_phy_kpi",100, 0.0, 3.0);
  TH1D *h_cl = new TH1D("h_cl","h_cl",100,0.0, 1.0);
  
  /////////////////////////////////////////////////////////////////////////////////////////
  //SANTORO PLOTS
  //
  TH2D *h2_kp_phitheta = new TH2D("h2_kp_phitheta","h2_kp_phitheta", 200, 0.0, 60.0, 200, -30.0, 30.0 );
  
  TH2D *h_mmekx_kmntm = new TH2D("h_mmekx_kmntm","h_mmekx_kmntm",300, 0.0, 12.6, 300, 1.2, 2.75);
  TH2D *h_mmekx_kmntm2 = new TH2D("h_mmekx_kmntm2","h_mmekx_kmntm2",300, 0.0, 12.6, 300, 1.2, 2.75);
  TH2D *h_imkk_kktheta = new TH2D("h_imkk_kktheta","h_imkk_kktheta",300, 0.0, 60.0, 300, 0.7, 2.0);
  TH2D *h_imkk_kktheta2 = new TH2D("h_imkk_kktheta2","h_imkk_kktheta2",300, 0.0, 60.0, 300, 0.7, 2.0);

  TH2D *h_prperp_mmekx = new TH2D("h_prperp_mmekx","h_prperp_mmekx",300,0.0,18, 300, 0.0, 60.0);
  TH2D *h_mmekx_test = new TH2D("h_mmekx_test","h_mmekx_test",300, 0,2, 300, -2.0, 2.0);
  TH2D *h_mm_pion_sep = new TH2D("h_mm_pion_sep","h_mm_pion_sep",300,0.0,1.4, 300, 0.0, 1.4);

  
  
  ////////////////////////////////////////////////////////////////////////////////////////
  //KINFIT PLOTS
  TH2D *h_e_p = new TH2D("h_e_p","h_e_p",100,0.0, 6.0, 200, 0.0, 60.0);
  TH2D *h_pr_p = new TH2D("h_pr_p","h_pr_p",100,0.0, 6.0, 200, 0.0, 60.0);
  TH2D *h_kp_p = new TH2D("h_kp_p","h_kp_p",100,0.0, 4.0, 200, 0.0, 60.0);
  TH2D *h_km_p = new TH2D("h_km_p","h_km_p",100,0.0, 4.0, 200, 0.0, 60.0);
  TH1D *h_phy_ekX_kf = new TH1D("h_phy_ekX_kf","h_phy_ekX_kf", 200, 0.9, 1.6);
  TH1D *h_phy_imkk_kf_b = new TH1D("h_phy_imkk_kf_b","h_phy_imkk_kf_b", 200, 0.8, 1.5);
  TH1D *h_phy_imkk_kf_a = new TH1D("h_phy_imkk_kf_a","h_phy_imkk_kf_a", 200, 0.8, 1.5);
    
  TLorentzVector lv_ebeam;
  lv_ebeam.SetPxPyPzE( 0, 0, 5.498, 5.498);
  TLorentzVector lv_target;
  lv_target.SetPxPyPzE( 0, 0, 0, 0.938);

  int hel = 0;

  //KinFitter *kinFit = nullptr;
  //TString fname   = "fittedEvents.root";
  //TFile *fileout  = new TFile( fname, "recreate" );
  //TTree *treeout = new TTree("kinfit", "Kin. Fitted Events");
  
  //TBranch* b_kinfit = treeout->Branch("KinFitEvents", "KinFitter", &kinFit);
  //TBranch* b_hel    = treeout->Branch("hel"   , &hel   );
  
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
    TLorentzVector lv_lambda = lv_pr + lv_km;
    
    TLorentzVector lv_ekX = lv_ebeam + lv_target - lv_el - lv_kp;

    TLorentzVector lv_epkX = lv_ebeam + lv_target - lv_el - lv_pr - lv_kp;

    TLorentzVector lv_mod_kp;
    lv_mod_kp.SetPxPyPzE( lv_kp.Px(), lv_kp.Py(), lv_kp.Pz(),  sqrt(lv_kp.Vect().Mag2() + (0.139*0.139)) );
    TLorentzVector mod_epkX = lv_ebeam + lv_target - lv_el - lv_pr - lv_mod_kp; 
    //h_mm_pion_sep->Fill( lv_epkX.M(), mod_epkX.M() );
    
    /////////////////////////////////////////////////////////
    //
    //std::cout << " >> " << lv_pr.M() << " " << lv_kp.M() << std::endl;
    if( topology != 0 ){
      TLorentzVector *P_e = ((ex_event->electron));
      TLorentzVector *P_p = ((ex_event->proton));
      TLorentzVector *P_kp = ((ex_event->kaon_pos));
      TLorentzVector *P_km = ((ex_event->kaon_neg));
      
      std::vector<Double_t> sigmas_e, sigmas_p, sigmas_kp, sigmas_km;
      //sigmas_e    = getResErr_DC( *P_e );
      //sigmas_p    = getResErr_DC( *P_p );
      //sigmas_kp = getResErr_DC( *P_kp );
      //sigmas_km = getResErr_DC( *P_km );
    
      //sigmas_p[0]*=5;
      //sigmas_p[1]*=20;
      //sigmas_p[2]*=10;

      
      //sigmas_kp[0]*=5;
      //sigmas_kp[1]*=20;
      //sigmas_kp[2]*=10;

      //sigmas_km[0]*=5;
      //sigmas_km[1]*=20;
      //sigmas_km[2]*=10;
      //*/    

      /*   TLorentzVector P_phi = (lv_ebeam) + (lv_target) - (*P_e + *P_p);
      h_phy_imkk_kf_b->Fill(P_phi.M());
      double phimass = P_phi.M();
      P_phi.SetVectM( P_phi.Vect(), 1.019);

      P_km->SetVectM( P_km->Vect(), 0.493);
      
      KFParticle  KF_e     =  KFParticle(*P_e, sigmas_e,"e");
      KFParticle  KF_p     =  KFParticle(*P_p, sigmas_p,"p");
      KFParticle  KF_kp     =  KFParticle(*P_kp, sigmas_km,"kp");
      KFParticle  KF_km     =  KFParticle(*P_km, {} ,"km", 1, 0 );
      
      //KFParticle  KF_phi    =  KFParticle(P_phi, {}, "phi", 1, 0);
      
      //KFParticle  KF_km    =  KFParticle(*P_km, sigmas_km,"km");

      std::vector<KFParticle> KFParticles;
      KFParticles.push_back(KF_e);
      KFParticles.push_back(KF_p);
      //KFParticles.push_back(KF_phi);
      KFParticles.push_back(KF_kp);
      KFParticles.push_back(KF_km);
	
      ////, KF_kp, KF_km };
      kinFit = new KinFitter( KFParticles ); // Initiali Kin. Fit with list of particles
      kinFit->SetM_targ(0.938);              // Set the target to be the proton (Default is helium, M_targ = 3.74 GeV/c^2 )
      kinFit->SetBeam(5.498);
      kinFit->DoFitting(100);                // Do the fit 100 times or until chi^2 increases

      double cl_kinfit = kinFit->confLevel;
      //std::cout << " >> " << cl_kinfit << std::endl;
	
      if ( cl_kinfit >= 0.3 ){
	TLorentzVector lvkf_e = kinFit->KFParticlesOut[0].P;
	TLorentzVector lvkf_pr = kinFit->KFParticlesOut[1].P;
	TLorentzVector lvkf_kp = kinFit->KFParticlesOut[2].P;
	//TLorentzVector lvkf_phi = kinFit->KFParticlesOut[2].P;
	//TLorentzVector lvkf_kp = kinFit->KFParticlesOut[2].P;
	TLorentzVector lvkf_km = kinFit->KFParticlesOut[3].P;
	TLorentzVector lvkf_phi = lvkf_kp + lvkf_km;
	
	
	h_e_p->Fill( lvkf_e.P(), lvkf_e.Theta()*180.0/3.14159265358 );
	h_pr_p->Fill( lvkf_pr.P(), lvkf_pr.Theta()*180.0/3.14159265358  );
	h_kp_p->Fill( lvkf_phi.P(), lvkf_phi.Theta()*180.0/3.14159265358  );
	//h_km_p->Fill( lvkf_km.P(), lvkf_km.Theta()*180.0/3.14159265358  );
	
	h_phy_imkk_kf_a->Fill( lvkf_phi.M() );
	//    TLorentzVector lvkp_phi =
	treeout->Fill();
	//}
      }
    
      h_cl->Fill(cl_kinfit);
      
      //}
      delete kinFit;
    
      //*/
      /////////////////////////////////////////////////////////    
    }
    TLorentzVector lv_kppi;
    lv_kppi.SetPxPyPzE( lv_kp.Px(), lv_kp.Py(), lv_kp.Pz(),  (0.139*0.139) );
    
    TLorentzVector lv_ekpiXmass  = lv_ebeam + lv_target - lv_el - lv_pr - lv_kppi;   
    double t_pr = (lv_pr - lv_target).M2();
    
    TLorentzVector lv_vph = lv_ebeam - lv_el;
    TLorentzVector lv_b_el(lv_el), lv_b_pr(lv_pr), lv_b_kp(lv_kp), lv_b_km(lv_km), lv_b_phi(lv_ebeam + lv_target - lv_el - lv_pr), lv_b_lambda(lv_ekX);
    
    
    TLorentzVector lv_wsys = -(lv_vph + lv_target);
    //TVector3 v_boost = lv_wsys.BoostVector();
    TVector3 v_boost = -lv_wsys.BoostVector();
    //lv_vph.Boost(v_boost);
    
    lv_b_el.Boost(v_boost);
    lv_b_pr.Boost(v_boost);
    lv_b_phi.Boost(v_boost);
    lv_b_lambda.Boost(v_boost);
    lv_b_kp.Boost(v_boost);
    lv_b_km.Boost(v_boost);
      
      
    double ang_kk = (180.0/3.14159265) * (lv_kp.Vect().Angle(lv_km.Vect()) );
    double ang_kpr = (180.0/3.14159265) * (lv_b_kp.Angle(lv_b_lambda.Vect()) );

                 
      
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


    h_phy_ekX->Fill( lv_ekX.M() );
    if( !(lv_epkX.M() <= 0.2 && mod_epkX.M() <= 0.4) ){
      h_phy_ekX_ac->Fill( lv_ekX.M() );
      h_mm_pion_sep->Fill( lv_epkX.M(), mod_epkX.M() );
    }


    
    if( topology == 0 ){
      TLorentzVector lv_phi = lv_kp + lv_km;

      
      h_phy_ang_kk->Fill(ang_kk);
      //h_phy_ekX->Fill( lv_ekX.M() );
      h_phy_epkXIM->Fill( lv_phi.M(), lv_km.M() );
      h_phy0_mm2->Fill(lv_km.M());

      
      
      if( (lv_ekX.M() <= 1.53 || lv_ekX.M() >= 1.57) && (lv_phi.M() <= 1.04 && lv_phi.M() >= 0.99) ){
	h_phy_kpi->Fill(lv_kppi.M());
	
	//if(
	/*	if(helicity > 0){
	  h_asy_phi_P->Fill(n_ang);
	}
	else{
	  h_asy_phi_N->Fill(n_ang);
	}
	*/
	
      }
      
    }
    if( topology == 0 ){
      TLorentzVector lv_phi = lv_kp + lv_km;
      //double ang_kk = (180.0/3.14159265) * (lv_kp.Vect().Angle(lv_km.Vect()) );    
      //std::cout << lv_b_kp.Theta()*(180.0/3.14159265) << std::endl;
      //std::cout << " >> KM " << lv_b_km.Theta()*(180.0/3.14159265) << std::endl;
      
      h_phy_ang_decay->Fill(ang_kk,ang_kpr);

      h_phy3_ang_kk->Fill( ang_kk );
      h_phy3_imkk->Fill(lv_phi.M());
      h_phy3_mm2->Fill(lv_km.M());
      h_phy_angphi->Fill(n_ang);

      h_phy03_imkk_kpr->Fill( lv_phi.M(), lv_ekX.M() );

      // if( helicity > 0 ){
      //	h_asy_phi_P->Fill(n_ang);
      //}
      //else if ( helicity < 0 ) {
      //	h_asy_phi_N->Fill(n_ang);
      //}
      /*if(helicity > 0){
	h_asy_phi_P->Fill(n_ang);
      }
      else{
	h_asy_phi_N->Fill(n_ang);
	} */     
    }
    TLorentzVector lv_phi_all = lv_kp + lv_km;
    double kk_angle = (180.0/3.14159265) * (lv_kp.Vect().Angle(lv_km.Vect()) );
    double kkang_cos = TMath::Sin((lv_kp.Vect().Angle(lv_km.Vect()) ));
    double phi_pr_angle = (180.0/3.14159265) *( lv_phi_all.Vect().Angle( lv_pr.Vect() ));
    double pr_perp = lv_pr.Perp();//Vect().Mag();
    h_prperp_mmekx->Fill(pr_perp, phi_pr_angle );
  
    if( lv_ekX.M() <= 1.50 && lv_ekX.M() >= 1.549 ){
      h_imkk_kktheta->Fill(kk_angle, lv_ekX.M());
      h_mmekx_kmntm->Fill(lv_km.E(), lv_ekX.M());
      h_phy03_imkk->Fill( lv_phi_all.M() );	  

    }
      h_imkk_kktheta2->Fill(kk_angle, lv_ekX.M());
      
    
    h_phy_t->Fill(-t_pr);
    h_mmekx_kmntm2->Fill( lv_km.E() , lv_ekX.M() );

    h2_kp_phitheta->Fill( (180.0/3.14159265) * lv_kp.Vect().Theta(), (180.0/3.14159265) * lv_kp.Vect().Phi() );
    
    double deltammepkX = lv_ekX.M() - 1.0194;
    //    if( deltammepkX <= 0.48 || deltammepkX >= 0.52 ){

    if( (lv_ekX.M() <= 1.51 || lv_ekX.M() >= 1.57) && (lv_phi_all.M() <= 1.04 && lv_phi_all.M() >= 0.99) ){
      
      h_mmekx_test->Fill( kkang_cos, deltammepkX);
      //h_phy03_imkk_kpr->Fill( lv_phi_all.M(), lv_ekX.M() );
    }

    //h_mm_pion_sep->Fill( lv_epkX.M() , lv_ekpiXmass.M());
    h_phy_kpi->Fill( lv_ekpiXmass.M() );

    h_phy03_imkk_kf->Fill( lv_phi_all.M() );
    //h_phy03_imkk->Fill( lv_phi_all.M() );
      
    double imkk_cut_top = 0.51303;
    double imkk_cut_bot = 0.472634;
    if( (lv_km.M() <= imkk_cut_top && lv_km.M() >= imkk_cut_bot)
 	&& !(lv_epkX.M() <= 0.2 && mod_epkX.M() <= 0.4) 
	&& (lv_ekX.M() <= 1.495 || lv_ekX.M() >= 1.535) 
	&& (lv_phi_all.M() >= 1.0 && lv_phi_all.M() <= 1.034) ){

      //h_phy03_imkk_kpr->Fill( lv_phi_all.M(), lv_ekX.M() );
      

      if(helicity > 0){
	h_asy_phi_P->Fill(n_ang);
      }
      else{
	h_asy_phi_N->Fill(n_ang);
      }
      
      
    }
    
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

  TH1D *h_asy_delta_top = new TH1D( "h_asy_delta_top", "h_asy_phi_top", 10, -180.0, 180.0);// = h_asy_phi_P;
  TH1D *h_asy_delta_bot = new TH1D( "h_asy_delta_bot", "h_asy_phi_bot", 10, -180.0, 180.0);;// = h_asy_phi_P;
  
  h_asy_delta_top->Add(h_asy_delta_P, h_asy_delta_N, -1.0);
  h_asy_delta_bot->Add(h_asy_delta_P, h_asy_delta_N, 1.0);
  
  TH1D *h_asy_delta_final =  new TH1D("h_asy_delta_final","h_asy_delta_final",10,-180.0, 180.0);
  h_asy_delta_final->Divide(h_asy_delta_top, h_asy_delta_bot, 1.0, 1.0);
  
  TGraph *g_asy_delta_final = new TGraph(h_asy_delta_final);
  g_asy_delta_final->SetName("g_asy_delta_final");
  g_asy_delta_final->SetTitle("#phi BSA");

  
  TCanvas *c_mm2 = new TCanvas("c_mm2","c_mm2",900,900);
  //c_mm2->Divide(1,1);
  //c_mm2->cd(1);
  h_phy_ekX->GetXaxis()->SetTitle("M.M. ek^{+}X [GeV]");
  h_phy_ekX->GetXaxis()->CenterTitle();
  //h_phy_ekX->SetFillStyle(3375);
  h_phy_ekX->SetFillColorAlpha(kRed-4, 0.35);
  h_phy_ekX->Draw();
  //h_phy_ekX_ac->SetFillStyle(3357);
  h_phy_ekX_ac->SetFillColorAlpha(kBlue-4, 0.35);
  h_phy_ekX_ac->Draw("same");
  //c_mm2->cd(2);
  //h_phy3_mm2->Draw();

  TCanvas *c2 = new TCanvas("c2","c2",1600,800);
  //gStyle->SetOptStat(00000);
  //gStyle->SetOptFit(0011111);
  c2->Divide(2,1);
  c2->cd(1);
  h_phy0_mm2->GetXaxis()->SetTitle("M.M. epk^{+}X [GeV]");
  h_phy0_mm2->GetXaxis()->CenterTitle();
  //TF1 *f_km_pol = new TF1("f_km_pol","[0] + [1]*x + [2]*x*x", 0.505, 0.69);
  //h_phy0_mm2->Fit("f_km_pol","O+R+");
  
  /*TF1 *f_km_mass = new TF1("f_km_mass","([0]*exp(-0.5*((x-[1])/[2])^2))", 0.45, 0.505);// + ([3] + [4]*x + [5]*x*x)", 0.45, 0.7 );//0.4701, 0.52);
  //f_km_mass->SetRange(0.45,0.54);
  f_km_mass->SetParameter(0,6500);
  f_km_mass->SetParameter(1,0.4976);
  f_km_mass->SetParameter(2,0.05);
  h_phy0_mm2->Fit("f_km_mass","O+R+");    
  */
		      
  TF1 *f_km_mass_tot = new TF1("f_km_mass_tot","([0]*exp(-0.5*((x-[1])/[2])^2)) + ([3] + [4]*x + [5]*x*x)", 0.45, 0.7);
  f_km_mass_tot->SetRange(0.45, 0.7);
  f_km_mass_tot->SetParameter(0, 6500);//f_km_mass->GetParameter(0));
  f_km_mass_tot->SetParameter(1, 0.4976);//f_km_mass->GetParameter(1));
  f_km_mass_tot->SetParameter(2, 0.05);//f_km_mass->GetParameter(2));  
  f_km_mass_tot->SetParameter(3, 1.0);//f_km_pol->GetParameter(0));
  f_km_mass_tot->SetParameter(4, 1.0);//f_km_pol->GetParameter(1));
  f_km_mass_tot->SetParameter(5, 1.0);//f_km_pol->GetParameter(2));
  h_phy0_mm2->Fit("f_km_mass_tot","R+");    
  
  h_phy0_mm2->Draw();
  f_km_mass_tot->Draw("same");
  //f_km_mass->Draw("same");
  //f_km_pol->Draw("same");
 
  c2->cd(2);

  double sig_fit = 1.5;
  double mm_mean = f_km_mass_tot->GetParameter(1);
  double mm_cut_top = f_km_mass_tot->GetParameter(1) + sig_fit*f_km_mass_tot->GetParameter(2);
  double mm_cut_bot = f_km_mass_tot->GetParameter(1) - sig_fit*f_km_mass_tot->GetParameter(2);
  std::cout << " >> " << mm_mean << " " << mm_cut_top << " " << mm_cut_bot << std::endl;
  
  TLine *cut_top = new TLine(0.8, mm_cut_top, 1.3, mm_cut_top);
  TLine *cut_mean = new TLine(0.8, mm_mean, 1.3, mm_mean);
  TLine *cut_bot = new TLine(0.8, mm_cut_bot, 1.3, mm_cut_bot);
  cut_top->SetLineColor(kRed);
  cut_mean->SetLineColor(kRed);
  cut_bot->SetLineColor(kRed);
  cut_top->SetLineWidth(2);
  cut_mean->SetLineWidth(2);
  cut_mean->SetLineStyle(2);
  cut_bot->SetLineWidth(2);
    
  h_phy_epkXIM->GetXaxis()->SetTitle("IM_{kk} [GeV]");
  h_phy_epkXIM->GetXaxis()->CenterTitle();
  h_phy_epkXIM->GetYaxis()->SetTitle("M.M. epk^{+}X [GeV]");
  h_phy_epkXIM->GetYaxis()->CenterTitle();
  h_phy_epkXIM->Draw("colz");
  cut_bot->Draw("same");
  cut_mean->Draw("same");
  cut_top->Draw("same");  
  
  TCanvas *c4 = new TCanvas("c4","c4",900,900);
  //gPad->SetLogz();
  TLine *kp_bot = new TLine(0.70, 1.495, 2.0 , 1.495);
  TLine *kp_top = new TLine(0.70, 1.535, 2.0 , 1.535);
  kp_bot->SetLineColor(kRed);
  kp_top->SetLineColor(kRed);
  kp_bot->SetLineWidth(2);
  kp_top->SetLineWidth(2);
  h_phy03_imkk_kpr->GetXaxis()->SetTitle("I.M_{kk} [GeV]");
  h_phy03_imkk_kpr->GetXaxis()->CenterTitle();
  h_phy03_imkk_kpr->GetYaxis()->SetTitle("M.M. ek^{+}X [GeV]");
  h_phy03_imkk_kpr->GetYaxis()->CenterTitle();
  h_phy03_imkk_kpr->Draw("colz");
  kp_bot->Draw("same");
  kp_top->Draw("same");
   
  
  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  TF1 *f_asy = new TF1("f_asy","[0]*sin(x*(3.1415926/180.0)) + [1]", -180.0, 180.0);
		      

  g_asy_phi_final->GetXaxis()->SetTitle("#phi [deg]");
  g_asy_phi_final->GetXaxis()->CenterTitle();
  g_asy_phi_final->GetYaxis()->SetTitle("BSA");
  g_asy_phi_final->GetYaxis()->CenterTitle();
  g_asy_phi_final->SetMarkerStyle(23);
  g_asy_phi_final->SetMarkerSize(1);
  g_asy_phi_final->Fit(f_asy);
  g_asy_phi_final->Draw("AP");
  f_asy->Draw("same");

  TCanvas *c5 = new TCanvas("c5","c5",900,900);
  //c5->Divide(2,1);
  //c5->cd(1);
  TF1 *f_phi_gaus = new TF1("f_phi_gaus","([0]*exp(-0.5*((x-[1])/[2])^2)) + ([3]*sqrt(x*x - 0.98*0.98) + [4]*(x*x - 0.98*0.98))", 0.98, 1.15);// + ([3]*sqrt(x*x - [4]) + [5](x*x - [4]))", 0.95, 1.15);
  TF1 *f_phi_pol = new TF1("f_phi_pol", "([0]*sqrt(x*x - 0.98*0.98) + [1]*(x*x - 0.98*0.98))", 0.98, 1.15);
  
  TF1 *f_phi_bw = new TF1("myfunc", "TMath::BreitWigner(x,1.019,0.00426)", 0.98, 1.15);
  f_phi_bw->SetLineColor(kBlue-2);
  f_phi_gaus->SetParameter(0,h_phy03_imkk->GetMaximum());
  f_phi_gaus->SetParameter(1,1.019);
  f_phi_gaus->SetParameter(2,0.0081);
  f_phi_gaus->SetParameter(3, 1.0);
  f_phi_gaus->SetParameter(4, 1.0);
  
  h_phy03_imkk->Fit(f_phi_gaus,"R+");
  h_phy03_imkk->SetMarkerStyle(28);
  h_phy03_imkk->SetMarkerSize(1.5);
  h_phy03_imkk->GetXaxis()->SetTitle("I.M. k^{+}k^{-} [GeV]");
  h_phy03_imkk->GetXaxis()->CenterTitle();
  h_phy03_imkk->Draw("E");

  f_phi_pol->SetParameter(0, f_phi_gaus->GetParameter(3));
  f_phi_pol->SetParameter(1, f_phi_gaus->GetParameter(4));
  f_phi_pol->SetLineStyle(2);
  
  f_phi_pol->Draw("same");
  f_phi_gaus->Draw("same");
  f_phi_bw->Draw("same");
  //h_mmekx_kmntm->Draw("colz");

  ////////////////////
  //YIELD CALCULATION
  TAxis *axis = h_phy03_imkk->GetXaxis();
  int bmin = axis->FindBin(1.0);
  int bmax = axis->FindBin(1.044);
  double signal = h_phy03_imkk->Integral(bmin,bmax);
  //signal -= h_phy03_imkk->GetBinContent(bmin)*(1.0 - axis->GetBinLowEdge(bmin))/(axis->GetBinWidth(bmin));
  //signal -= h_phy03_imkk->GetBinContent(bmax)*(axis->GetBinUpEdge(bmin) - 1.044 )/(axis->GetBinWidth(bmax));
  std::cout << " signal " << " 1.0 - " << axis->GetBinLowEdge(bmin) << std::endl;
  std::cout << " signal " << " -1.045 + " << axis->GetBinUpEdge(bmin) << std::endl;

  std::cout<< " >> SIGNAL " << signal << std::endl;
  
  double bin_size = (h_phy03_imkk->GetXaxis()->GetXmax() - h_phy03_imkk->GetXaxis()->GetXmin() )/((double)(h_phy03_imkk->GetXaxis()->GetNbins()));
  //std::cout << bin_size << std::endl;
  double phi_yield = (sqrt(2.0 * TMath::Pi())/bin_size) * (f_phi_gaus->GetParameter(0)*f_phi_gaus->GetParameter(2));///h_phy03_imkk.;
  std::cout << " PHI YIELD IS USING SANTORO METHOD IS " << phi_yield << std::endl;
  std::cout << " INTEGRATING BACKGROUND BETWEEN BIN " << bmin << " -> " << bmax << std::endl;
  double tot_background = 0.0;
  double b0 = f_phi_gaus->GetParameter(3);
  double b1 = f_phi_gaus->GetParameter(4);
  
  for( int b = bmin; b <= bmax; b++ ){
    double x = axis->GetBinCenter(b);
    double bg = (b0*sqrt(x*x - 0.98*0.98) + b1*(x*x - 0.98*0.98));// * axis->GetBinWidth(b);
    std::cout << " x " << x << " bg " << bg << " " << axis->GetBinWidth(b) << std::endl;
    tot_background+=bg;
  }
  std::cout << " TOTAL BACKGROUND IS " << tot_background << std::endl;
  std::cout << " FINAL PHI RATE IS " << signal - tot_background << std::endl;
  
  //c5->cd(2);
  //h_phy03_imkk_kf->Draw();
  
  //h_mmekx_kmntm2->Draw("colz");

  TCanvas *c6 = new TCanvas("c6","c6",900,900);
  c6->Divide(2,3);
  c6->cd(1);
  h_e_p->Draw("colz");
  c6->cd(2);
  h_pr_p->Draw("colz");
  c6->cd(3);
  h_kp_p->Draw("colz");
  c6->cd(4);
  h_km_p->Draw("colz");
  c6->cd(5);
  gPad->SetLogy();
  h_cl->Draw();
 
  TCanvas *c7 = new TCanvas("c7","c7",900,900);
  c7->Divide(2,2);
  c7->cd(1);
  h_phy_ekX_kf->Draw();
  c7->cd(2);
  TLine *l_phimass_b = new TLine( 1.019, 0.0, 1.019, h_phy_imkk_kf_b->GetMaximum() );
  h_phy_imkk_kf_b->Draw();
  l_phimass_b->Draw("same");
  c7->cd(3);
  TLine *l_phimass_a = new TLine( 1.019, 0.0, 1.019, h_phy_imkk_kf_b->GetMaximum() );
  h_phy_imkk_kf_a->Draw();
  l_phimass_a->Draw("same");
  c7->cd(4);
  h_phy_kpi->Draw();


  TCanvas *c8 = new TCanvas("c8","c8",900,900);
  //gStyle->SetOptStat(00000);
  TLine *pion_cut_x = new TLine(0.2, 0.0, 0.2, 1.4);
  TLine *pion_cut_y = new TLine(0.0, 0.4, 1.4, 0.4);
  h_mm_pion_sep->GetXaxis()->SetTitle("M.M. epk^{+}X [GeV]");
  h_mm_pion_sep->GetYaxis()->SetTitle("M.M. epk^{+}_{#pi mass}X [GeV]");
  h_mm_pion_sep->GetXaxis()->CenterTitle();
  h_mm_pion_sep->GetYaxis()->CenterTitle();
  h_mm_pion_sep->Draw("colz");
  pion_cut_x->SetLineWidth(2);
  pion_cut_y->SetLineWidth(2);
  pion_cut_x->SetLineColor(kRed);
  pion_cut_y->SetLineColor(kRed);
  pion_cut_x->Draw("same");
  pion_cut_y->Draw("same");

  //fileout->Write();
  //fileout->Close();

  TCanvas *c9 = new TCanvas("c9","c9",900,900);
  h_phy_ang_decay->Draw("colz");
  
}



//double getPhiAngle( TLorentzVector temp_el, TLorentzVector temp_pr, TLorentzVector
