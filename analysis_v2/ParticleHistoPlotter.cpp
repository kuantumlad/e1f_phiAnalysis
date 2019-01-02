#include "ParticleHistoPlotter.h"


#include "TH1D.h"
#include "TH2D.h"

ParticleHistoPlotter::ParticleHistoPlotter(){

    
  
}

ParticleHistoPlotter::~ParticleHistoPlotter(){

  

}


void ParticleHistoPlotter::initializeHistos(const char* temp_part_name ){

  f  = new TFile(Form("h1_%s.root",temp_part_name),"RECREATE");

  
  h_p = new TH1D(Form("h_p_%s",temp_part_name),Form("h_p_%s",temp_part_name), 100, 0.0, 0.0);
  h_theta = new TH1D(Form("h_theta_%s",temp_part_name),Form("h_theta_%s",temp_part_name), 100, 0.0, 0.0);
  h_phi = new TH1D(Form("h_phi_%s",temp_part_name),Form("h_phi_%s",temp_part_name), 300, 0.0, 0.0);

  h2_ptheta = new TH2D(Form("h2_ptheta_%s",temp_part_name),Form("h_ptheta_%s",temp_part_name), 100, 0.0, 5.9, 100, 0.0, 0.0);
  h2_pphi = new TH2D(Form("h2_pphi_%s",temp_part_name),Form("h_pphi_%s",temp_part_name), 100, 0.0, 5.9, 300, 0.0, 0.0);
  h2_thetaphi = new TH2D(Form("h2_thetaphi_%s",temp_part_name),Form("h_thetaphi_%s",temp_part_name), 100, 0.0, 0.0, 300, 0.0, 0.0);

  
}


void ParticleHistoPlotter::RescaleAxis(){

  
   
}

void ParticleHistoPlotter::DisplayHistos( const char* temp_part_name ){

  
  TCanvas *c1 = new TCanvas(Form("h1_%s",temp_part_name),Form("h1_%s",temp_part_name), 900,900);
  h_p->SetTitle("Momentum");
  h_p->GetXaxis()->SetTitle("momentum [GeV]");
  h_p->Draw();
  c1->Print(Form("h1_%s.pdf(",temp_part_name),"pdf");
  h_theta->Draw();
  c1->Print(Form("h1_%s.pdf",temp_part_name),"pdf");
  h_phi->Draw();
  c1->Print(Form("h1_%s.pdf)",temp_part_name),"pdf");

  TCanvas *c2 = new TCanvas(Form("h2_%s",temp_part_name),Form("h2_%s",temp_part_name), 900,900);
  gStyle->SetOptStat(000000);
  h2_ptheta->SetTitle("Momentum vs Theta");
  h2_ptheta->GetXaxis()->SetTitle("Momentum [GeV]");
  h2_ptheta->GetXaxis()->CenterTitle();
  h2_ptheta->GetYaxis()->SetTitle("#theta [deg]");  
  h2_ptheta->GetYaxis()->CenterTitle();
  h2_ptheta->Draw("colz");
  c2->Print(Form("h2_%s.pdf(",temp_part_name),"pdf");
  h2_pphi->SetTitle(" Momentum vs #phi");
  h2_pphi->GetXaxis()->SetTitle("Momentum [GeV]");
  h2_pphi->GetXaxis()->CenterTitle();
  h2_pphi->GetYaxis()->SetTitle("#theta [deg]");
  h2_pphi->GetYaxis()->CenterTitle();
  h2_pphi->Draw("colz");
  c2->Print(Form("h2_%s.pdf",temp_part_name),"pdf");
  h2_thetaphi->SetTitle("#theta vs #phi");
  h2_thetaphi->GetXaxis()->SetTitle("#theta [deg]");
  h2_thetaphi->GetXaxis()->CenterTitle();
  h2_thetaphi->GetYaxis()->SetTitle("#phi [deg]");
  h2_thetaphi->GetYaxis()->CenterTitle();   
  h2_thetaphi->Draw("colz");
  c2->Print(Form("h2_%s.pdf)",temp_part_name),"pdf");

  f->Write();
  f->Save();
  f->Close();
  
}




