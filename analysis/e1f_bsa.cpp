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
#include "GetBSA.h"
#include "GetPhiAngle.h"
#include "GetBSAError.h"
#include "GetCutLimits.h"
#include "SetCutLimits.h"
#include "CalcAsyError.h"


TLorentzVector lv_ebeam;
TLorentzVector lv_target;

void e1f_bsa( ){
  
  TChain *fchain = new TChain("events");
  fchain->Add("/Users/bclary/Documents/work/thesis_work/e1f_phiAnalysis/clas-phi/data/events.root");
  ExclusiveEvents *ex_event = new ExclusiveEvents(fchain);
  Long64_t nentries = fchain->GetEntries();

  lv_ebeam.SetPxPyPzE( 0, 0, 5.498, 5.498);
  lv_target.SetPxPyPzE( 0, 0, 0, 0.938);

  double cut_lmb_max = 1.535;
  double cut_lmb_min = 1.495;
  double cut_phi_max = 1.034;
  double cut_phi_min = 1.0;

  double pr_cl_min = 0.1;//0.04;
  double kp_cl_min = 0.1;//0.05;
  
  int n_bins = 10;
  
  std::map< std::string, std::vector< std::vector<double> > > m_cutlimits; //element [0] = loose , element[1] = tight
  GetCutLimits(m_cutlimits);
  int n_cutrange = m_cutlimits["dist_cc"][0].size(); // same for all of the cuts
  
  std::vector< TH1D* > h_asy_pn;
  std::vector< TGraphErrors* > g_asy;
  std::vector<TF1* > f_asy_cuts;
  std::map< int, std::vector<TF1*> > m_f_asy_cuts;
  
  std::map< int, std::map<std::string, std::vector< double> > > cut_map_nom;
  std::map< int, std::map<std::string, std::vector< double> > > cut_map_1;
  std::map< int, std::map<std::string, std::vector< double> > > cut_map_2;
  SetCutLimits(0, cut_map_nom); // nom
  SetCutLimits(1, cut_map_1); // loose 
  SetCutLimits(2, cut_map_2); // tight 

  std::vector< std::map< int, std::map<std::string, std::vector< double > > > > v_cut_levels;
  v_cut_levels.push_back(cut_map_nom);
  v_cut_levels.push_back(cut_map_1);  
  v_cut_levels.push_back(cut_map_2);

  /*
    for( int p = 0; p < 3; p++ ){
    std::cout << " LEVEL " << p << std::endl;
    for( int j = 0; j < 9; j++ ){      
      std::cout << ">> CUT " << j << " MIN " << v_cut_levels[p][j]["dist_cc"][0] << " MAX " << v_cut_levels[p][j]["dist_cc"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_dcr1"][0] << " MAX " << v_cut_levels[p][j]["dist_dcr1"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_dcr3"][0] << " MAX " << v_cut_levels[p][j]["dist_dcr3"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_ecsf"][0] << " MAX " << v_cut_levels[p][j]["dist_ecsf"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_edep"][0] << " MAX " << v_cut_levels[p][j]["dist_edep"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_ecu"][0] << " MAX " << v_cut_levels[p][j]["dist_ecu"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_ecv"][0] << " MAX " << v_cut_levels[p][j]["dist_ecv"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_ecw"][0] << " MAX " << v_cut_levels[p][j]["dist_ecw"][1]  << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_vz"][0] << " MAX " << v_cut_levels[p][j]["dist_vz"][1] << std::endl;           
    }
    }
  */

  std::vector< std::string > v_cut_ntl = {"nom","loose","tight"};
  std::map<int, std::vector<TGraphErrors*> > m_g_asy;  
  for( int cut_tl = 0; cut_tl < 3; cut_tl++ ){
      std::cout<< " CUT STRENGTH: " <<  v_cut_ntl[cut_tl] << std::endl;      
    for( int cutlvl = 0; cutlvl < 9; cutlvl++ ){
      int p = cut_tl;
      int j = cutlvl;
      std::cout << ">> CUT " << j << " MIN " << v_cut_levels[p][j]["dist_cc"][0] << " MAX " << v_cut_levels[p][j]["dist_cc"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_dcr1"][0] << " MAX " << v_cut_levels[p][j]["dist_dcr1"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_dcr3"][0] << " MAX " << v_cut_levels[p][j]["dist_dcr3"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_ecsf"][0] << " MAX " << v_cut_levels[p][j]["dist_ecsf"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_edep"][0] << " MAX " << v_cut_levels[p][j]["dist_edep"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_ecu"][0] << " MAX " << v_cut_levels[p][j]["dist_ecu"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_ecv"][0] << " MAX " << v_cut_levels[p][j]["dist_ecv"][1] << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_ecw"][0] << " MAX " << v_cut_levels[p][j]["dist_ecw"][1]  << std::endl;
      std::cout << " CUT " << j << " MIN " << v_cut_levels[p][j]["dist_vz"][0] << " MAX " << v_cut_levels[p][j]["dist_vz"][1] << std::endl;

      
      h_asy_pn.push_back(new TH1D(Form("h_asy_P_%d_cutlvl_%d", cut_tl, cutlvl), Form("h_asy_P_%d_cutlvl_%d", cut_tl, cutlvl), n_bins, -180.0, 180.0 ));
      h_asy_pn.push_back(new TH1D(Form("h_asy_N_%d_cutlvl_%d", cut_tl, cutlvl), Form("h_asy_N_%d_cutlvl_%d", cut_tl, cutlvl), n_bins, -180.0, 180.0 ));
      f_asy_cuts.push_back(new TF1(Form("f_asy_%d_cutlvl_%d", cut_tl, cutlvl), "[0]*sin(x*(3.1415926/180.0))", -180.0, 180.0));
      
      for( Long64_t i = 0; i < nentries; i++ ){
	fchain->GetEntry(i);
	
	int helicity = ((ex_event->helicity));
	int topology = ((ex_event->topology));
      
	double pr_cl = (ex_event->alpha_kaon_pos);
	double kp_cl = (ex_event->alpha_proton);
      
	double dist_ecsf = (ex_event->dist_ecsf);
	double dist_ec_edep = (ex_event->dist_ec_edep);
	double dist_vz = (ex_event->dist_vz);
	double dist_cc_theta = (ex_event->dist_cc_theta);
	double dist_dcr1 = (ex_event->dist_dcr1);
	double dist_dcr3 = (ex_event->dist_dcr3);
	double dist_ecu = (ex_event->dist_ecu);
	double dist_ecv = (ex_event->dist_ecv);
	double dist_ecw = (ex_event->dist_ecw);
	double dist_cc = (ex_event->dist_cc);
	double imkk_cut_top = 0.51303;
	double imkk_cut_bot = 0.472634;

            
	if( ( v_cut_levels[cut_tl][cutlvl]["dist_cc"][0] <= dist_cc && dist_cc <= v_cut_levels[cut_tl][cutlvl]["dist_cc"][1]  ) &&
	    ( v_cut_levels[cut_tl][cutlvl]["dist_dcr1"][0] <= dist_dcr1 && dist_dcr1 <= v_cut_levels[cut_tl][cutlvl]["dist_dcr1"][1] ) &&
	    ( v_cut_levels[cut_tl][cutlvl]["dist_dcr3"][0] <= dist_dcr3 && dist_dcr3 <= v_cut_levels[cut_tl][cutlvl]["dist_dcr3"][1] ) &&
	    ( v_cut_levels[cut_tl][cutlvl]["dist_ecsf"][0] <= dist_ecsf && dist_ecsf <= v_cut_levels[cut_tl][cutlvl]["dist_ecsf"][1] ) &&
	    ( v_cut_levels[cut_tl][cutlvl]["dist_edep"][0] <= dist_ec_edep && dist_ec_edep <= v_cut_levels[cut_tl][cutlvl]["dist_edep"][1] ) &&
	    ( v_cut_levels[cut_tl][cutlvl]["dist_ecu"][0] <= dist_ecu && dist_ecu <= v_cut_levels[cut_tl][cutlvl]["dist_ecu"][1] ) &&
	    ( v_cut_levels[cut_tl][cutlvl]["dist_ecv"][0] <= dist_ecv && dist_ecv <= v_cut_levels[cut_tl][cutlvl]["dist_ecv"][1] ) &&
	    ( v_cut_levels[cut_tl][cutlvl]["dist_ecw"][0] <= dist_ecw && dist_ecw <= v_cut_levels[cut_tl][cutlvl]["dist_ecw"][1] ) &&
	    ( v_cut_levels[cut_tl][cutlvl]["dist_vz"][0] <= dist_vz && dist_vz <= v_cut_levels[cut_tl][cutlvl]["dist_vz"][1] ) ){	  
	
	  TLorentzVector lv_el = (*(ex_event->electron));       
	  TLorentzVector lv_pr = (*(ex_event->proton));
	  TLorentzVector lv_kp = (*(ex_event->kaon_pos));    
	  TLorentzVector lv_km = (*(ex_event->kaon_neg));
	  TLorentzVector lv_phi = lv_kp + lv_km;
	  TLorentzVector lv_kpr = lv_pr + lv_kp;
	  TLorentzVector lv_ekX = lv_ebeam + lv_target - lv_el - lv_kp;      
	  TLorentzVector lv_epkX = lv_ebeam + lv_target - lv_el - lv_pr - lv_kp;
	  TLorentzVector lv_mod_kp;
	  lv_mod_kp.SetPxPyPzE( lv_kp.Px(), lv_kp.Py(), lv_kp.Pz(),  sqrt(lv_kp.Vect().Mag2() + (0.139*0.139)) );
	  TLorentzVector mod_epkX = lv_ebeam + lv_target - lv_el - lv_pr - lv_mod_kp; 
	  
	  //--------------------------------------------------------------------
	  //IMPLEMENT BASIC MM CUT
	  if( pr_cl >= pr_cl_min && kp_cl >= kp_cl_min ){
	    if( (lv_ekX.M() <= cut_lmb_min || lv_ekX.M() >= cut_lmb_max)
		&& (lv_phi.M() <= cut_phi_max && lv_phi.M() >= cut_phi_min)
		&& (lv_km.M() <= imkk_cut_top && lv_km.M() >= imkk_cut_bot)
		&& !(lv_epkX.M() <= 0.2 && mod_epkX.M() <= 0.4) ){
	    
	      double phi_ang = GetPhiAngle(lv_ebeam, lv_target, lv_el, lv_pr, lv_kp, lv_km, lv_phi);
	    
	      //--------------------------------------------------------------------
	      //IMPLEMENT HELICITY SELECTION
	      if( helicity > 0 ){
		//h_asy_P->Fill(phi_ang);
		h_asy_pn[0]->Fill(phi_ang);
	      }
	      else{
		//h_asy_N->Fill(phi_ang);
		h_asy_pn[1]->Fill(phi_ang);
	      }       		
	    }
	  }
	}
      }
  
      double avg_pol = 0.75;	
      g_asy.push_back(GetBSA( n_bins, h_asy_pn[0], h_asy_pn[1], avg_pol ));
      h_asy_pn.clear();

      m_g_asy[cutlvl].push_back(GetBSA(n_bins, h_asy_pn[0], h_asy_pn[1], avg_pol ));
      
    }    
  }

  
  /*TCanvas *c3 = new TCanvas("c3","c3",800,800);
    TF1 *f_asy = new TF1("f_asy","[0]*sin(x*(3.1415926/180.0))", -180.0, 180.0);        
    g_bsa->GetXaxis()->SetTitle("#phi [deg]");
    g_bsa->GetXaxis()->CenterTitle();
    g_bsa->GetYaxis()->SetTitle("BSA");
    g_bsa->GetYaxis()->CenterTitle();
    g_bsa->SetMarkerStyle(23);
    g_bsa->SetMarkerSize(1);
    g_bsa->Fit(f_asy);
    g_bsa->Draw("AP");
    f_asy->Draw("same");
  */

  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  //TF1 *f_asy = new TF1("f_asy","[0]*sin(x*(3.1415926/180.0))", -180.0, 180.0);
  for( int c = 0; c < g_asy.size(); c++ ){
    g_asy[c]->GetXaxis()->SetTitle("#phi [deg]");
    g_asy[c]->GetXaxis()->CenterTitle();
    g_asy[c]->GetYaxis()->SetTitle("BSA");
    g_asy[c]->GetYaxis()->CenterTitle();
    g_asy[c]->SetMarkerStyle(24);
    g_asy[c]->SetMarkerSize(1);
    g_asy[c]->Fit(f_asy_cuts[c]);
    g_asy[c]->Draw("AP+same");
    f_asy_cuts[c]->SetLineColor(kViolet + c);
    f_asy_cuts[c]->Draw("same");
    std::cout << " >> FIT PARAMETER " << f_asy_cuts[c]->GetParameter(0) << std::endl;
  }

  map<int, std::vector<TGraphErrors*> >::iterator it;
  std::vector<std::string> cut_names;
  cut_names = {"dist_cc","dist_dcr1","dist_dcr3","dist_ecsf","dist_edep","dist_ecu","dist_ecv","dist_ecw","dist_vz"};    
  int c = 0;
  for( it = m_g_asy.begin(); it != m_g_asy.end(); ++it ){
    std::vector<TGraphErrors*>::iterator it2;
    std::cout << " CUT LVL " << it->first << std::endl;
    TCanvas *c4 = new TCanvas(Form("c%d",it->first),Form("c%d",it->first),800,800);
    std::vector<double> v_asy_cut_err;
    int d = 0;
    for( it2 = (it->second).begin(); it2 != (it->second).end(); ++it2 ){
      (*it2)->SetTitle(Form("BSA cut %s", cut_names[it->first].c_str()) );
      (*it2)->GetXaxis()->SetTitle("#phi [deg]");
      (*it2)->GetYaxis()->SetTitle("BSA");
      (*it2)->SetMarkerStyle(24);
      (*it2)->SetMarkerSize(1);
      (*it2)->Fit(f_asy_cuts[c]);
      (*it2)->Draw("AP+same");
      f_asy_cuts[c]->SetLineColor(kViolet + d);
      f_asy_cuts[c]->Draw("same");      
      double asy_cut = f_asy_cuts[c]->GetParameter(0);
      std::cout << " ASY " << asy_cut << std::endl;
      v_asy_cut_err.push_back(asy_cut);

      d++;     
      c++;      
    }
    std::cout << " err asy size " << v_asy_cut_err.size() << std::endl;
    std::cout << " CUT " << cut_names[it->first] << " ERROR: " <<  CalcAsyError( v_asy_cut_err ) << std::endl;
    v_asy_cut_err.clear();
  }
  
  //TCanvas *c4 = new TCanvas("c4","c4",800,800);
  //h_asy_P->SetLineColor(kRed);
  //h_asy_N->SetLineColor(kBlue);
  //h_asy_P->Draw();
  //h_asy_N->Draw("same");

  /*TCanvas *c5 = new TCanvas("c5","c5",800,800);
    TH1D *h_top = new TH1D("h_top","h_top",n_bins, -180.0, 180.0 );
    TH1D *h_bot = new TH1D("h_bot","h_bot",n_bins, -180.0, 180.0 );
    
    h_top->Add(h_asy_P, h_asy_N, 1.0, -1.0 );
    h_bot->Add(h_asy_P, h_asy_N, 1.0, 1.0 );
    h_top->SetLineColor(kRed);
    h_top->Draw();
    
    TCanvas *c6 = new TCanvas("c6","c6",800,800);
    h_bot->SetLineColor(kBlue);
    h_bot->Draw();
  */
    
    
}


