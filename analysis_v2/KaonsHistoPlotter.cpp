#include "KaonsHistoPlotter.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"

KaonsHistoPlotter::KaonsHistoPlotter(){



}


KaonsHistoPlotter::~KaonsHistoPlotter(){



}


void KaonsHistoPlotter::InitKaonHistos( const char* temp_part_name ){

  f = new TFile(Form("h_kaonpair_%s.root", temp_part_name ), "RECREATE");

  h_decayang = new TH1D(Form("h_decayang_%s",temp_part_name),Form("h_decayang_%s",temp_part_name), 100, 0.0, 65.0 );
  h_decayangp = new TH2D(Form("h_decayangp_%s",temp_part_name),Form("h_decayangp_%s",temp_part_name), 200, 0.0, 5.0, 200, 0.0, 65.0 );
  h_decayangp_prkm = new TH2D(Form("h_decayangp_prkm_%s",temp_part_name),Form("h_decayangp_prkm_%s",temp_part_name), 200, 0.0, 5.0, 200, 0.0, 65.0 );
  h_decayang_prkm_kpkm = new TH2D(Form("h_decayang_prkm_kpkm_%s",temp_part_name),Form("h_decayang_prkm_kpkm_%s",temp_part_name), 200, 0.0, 65.0, 200, 0.0, 65.0 );
  
}

void KaonsHistoPlotter::DisplayHistos( const char* temp_par ){

  TCanvas *c1 = new TCanvas(Form("h_kaons_%s",temp_par),Form("h_kaons_%s",temp_par), 900, 900 );
  gStyle->SetOptStat(00000);
  h_decayang->SetTitle("Decay Angle Between Kaons");
  h_decayang->GetXaxis()->SetTitle("Decay Angle in Lab Coord. [deg]");
  h_decayang->Draw("colz");
  c1->Print(Form("h_kaons_%s.pdf(",temp_par),"pdf");
  h_decayangp->SetTitle("Decay Angle of K^{-} + K^{+} vs Momentum");
  h_decayangp->GetXaxis()->SetTitle("Momentum [GeV]");
  h_decayangp->GetXaxis()->CenterTitle();
  h_decayangp->GetYaxis()->SetTitle("#theta [deg]");
  h_decayangp->GetYaxis()->CenterTitle();
  h_decayangp->Draw("colz");
  c1->Print(Form("h_kaons_%s.pdf",temp_par),"pdf");
  h_decayangp_prkm->SetTitle("Decay Angle of Pr + K^{+} vs P");
  h_decayangp_prkm->GetXaxis()->SetTitle("Momentum [GeV]");
  h_decayangp_prkm->GetXaxis()->CenterTitle();
  h_decayangp_prkm->GetYaxis()->SetTitle("Opening #theta [deg]");
  h_decayangp_prkm->GetYaxis()->CenterTitle();
  h_decayangp_prkm->Draw("colz");
  c1->Print(Form("h_kaons_%s.pdf",temp_par),"pdf"); 
  h_decayang_prkm_kpkm->SetTitle("Opening #theta_{Pr + K^{+}} vs Opening #theta_{K^{-} + K^{+}}");
  h_decayang_prkm_kpkm->GetXaxis()->SetTitle("Opening #theta_{K^{-} + K^{+}} [deg]");
  h_decayang_prkm_kpkm->GetXaxis()->CenterTitle();
  h_decayang_prkm_kpkm->GetYaxis()->SetTitle("Opening #theta_{Pr + K^{+}} [deg]");
  h_decayang_prkm_kpkm->GetYaxis()->CenterTitle();
  h_decayang_prkm_kpkm->Draw("colz");
  c1->Print(Form("h_kaons_%s.pdf)",temp_par),"pdf");
 
  
  f->Write();
  f->Save();
  f->Close();

}

