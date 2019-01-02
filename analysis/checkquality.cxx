#include "include/KinFitter.h"
#include "include/checkquality.h"
#include "include/physx/physx.h"

void checkquality(Double_t confCut = 0.15, Int_t nsig = 10000, Int_t nbg = 10000){
  // confCut : Where do you want to cut on the confidence level?
  
  // IO
  TString fname = "out/fittedEvents_" + to_string(nsig) + "_" + to_string(nbg) + ".root";
  TFile *filein = new TFile( fname , "read" );
  TTree *kintree = (TTree*)gROOT->FindObject("kinfit");
  
  fname = "out/histos_" + to_string(nsig) + "_" + to_string(nbg) + ".root"; 
  TFile *fileout = new TFile( fname, "recreate" );
  TTree * treeout = new TTree("histos", "Histograms for Kin. Fit");

  KinFitter * kinFit = nullptr;
  Int_t hel = 0;
  kintree->SetBranchAddress("KinFitEvents",  &kinFit);
  kintree->SetBranchAddress("hel",  &hel);

  // Set up histograms
  TH1D * h_confLevs_0 = new TH1D( "confLevs_0", "Confidence Levels; Conf. Lev. [ ]", 100, 0, 1 );
  TH1D * h_confLevs   = new TH1D( "confLevs", "Confidence Levels; Conf. Lev. [ ]", 100, 0, 1 );
  h_confLevs->SetFillColorAlpha( kBlue, 0.35 );

  TH1D * h_pulls    = new TH1D( "pulls", "Pull for ", 100, -5, 5 );
  std::vector<TH1D*> hh_pulls_0;
  std::vector<TH1D*> hh_pulls;

  //// Get the number of measured variables to set up pull histograms
  kintree->GetEntry(0);
  Int_t nvars_y = kinFit->nvars_y;
  for( Int_t ii = 0; ii < nvars_y; ii++ ){
    // Pull Distributions Before Confidence Level Cut (CLC)
    TH1D* hist = (TH1D*) h_pulls->Clone();  
    TString name = (TString) hist->GetName();
    TString title = (TString) hist->GetTitle();
    hist->SetName( name + "_" + to_string(ii) + "_0" );
    hist->SetTitle( title + " " + kinFit->names_y[ii] + "; " + kinFit->names_y[ii] + " " + kinFit->units_y[ii] );

    hh_pulls_0.push_back(hist);
    
    // Pull Distributions After CLC
    hist = (TH1D*) h_pulls->Clone();  
    name = (TString) hist->GetName();
    title = (TString) hist->GetTitle();
    hist->SetName( name + "_" + to_string(ii) );
    hist->SetTitle( title + "  " + kinFit->names_y[ii] + " (After CLC); " + kinFit->names_y[ii] + " " + kinFit->units_y[ii] );
    hist->SetFillColorAlpha( kBlue, 0.35 );

    hh_pulls.push_back(hist);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // Loop over events
  Int_t nevents = kintree->GetEntries();
  for( Int_t ievent = 0; ievent < nevents; ievent++){
    kintree->GetEntry(ievent);
    
    Double_t confLev = kinFit->confLevel;
    h_confLevs_0->Fill(confLev);


    ////////////////////////////////////
    //DO THIS TO GET PARTICLE LV OUT
    P_e = kinFit->KFParticlesOut[0].P;
    P_p = kinFit->KFParticlesOut[1].P;
    
    for( Int_t ii = 0; ii < nvars_y; ii++ ){
      Double_t pull = kinFit->pulls[ii];
      hh_pulls_0[ii]->Fill(pull);
    }
    // Fill if passes confidence level cut
    if( confLev > confCut ){

      
      h_confLevs->Fill(confLev);
      for( Int_t ii = 0; ii < nvars_y; ii++ ){
        Double_t pull = kinFit->pulls[ii];
        hh_pulls[ii]->Fill(pull);
      }
    }
      

    treeout->Fill();
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  //// Fit confidence level plateau to estimate signal to noise and change title with these stats
  TF1 * lfit = new TF1( "lfit", "pol0", 0.0, 1);
  h_confLevs->Fit("lfit", "", "", (1-confCut)/2, 1.0);
  Double_t plateau = lfit->GetParameter(0);
  Double_t snr    = calcSNR( h_confLevs, confCut, plateau );
  Double_t sigpct = calcSigPct( h_confLevs, confCut, plateau );
  
  stringstream ss;
  ss << fixed << setprecision(2) << confCut;
  clcut = ss.str();
  h_confLevs_0->SetTitle( ( string(h_confLevs_0->GetTitle()) + " ( " + to_string((Int_t)h_confLevs->GetEntries()) + " events, est. SNR = " + to_string(snr) + ", sig. pct. = " + to_string(sigpct) + "% with conf. cut @ " + clcut + ") ").c_str() );

    

  // Draw Histograms
  string imgoutdir = "images/";
  string pngName;
  TCanvas* c1 = new TCanvas;
  
  //// Draw Confidence Level Distribution
  pngName = "confLevs";
  gStyle->SetOptFit();
  h_confLevs_0->Draw();
  h_confLevs->Draw("same");
  c1->Update();
  ////// Draw horizontal and vertical lines representing plateau fit and CLC, respectively
  drawHLine( plateau, 1, 0.0, kRed, 1, 1); 
  drawVLine( confCut, 1.8*gPad->GetUymax()); 
  c1->SetLogy();
  printPNGandPDF( c1,  pngName, imgoutdir );

  //// Draw Pull Distributions
  pngName = "pulls";
  c1 = new TCanvas("c1","c1",3*400,nvars_y/3*400);
  c1->Divide(3 , nvars_y/3);
  for( Int_t ii = 0; ii < nvars_y; ii++ ){
    c1->cd(ii+1);
    hh_pulls[ii]->Draw();
    hh_pulls_0[ii]->Draw("same");
    
    // Actual fit to the data
    TF1 *gausFit2 = new TF1("gausFitt2", gausFitFunc, -5, 5, 3);
    gausFit2->SetParameters( hh_pulls[ii]->GetMaximum(), 0, 1);
    gausFit2->SetLineColorAlpha( kBlue, 0.35);
    gausFit2->SetParNames( "N (normalization)", "#mu (mean)", "#sigma (width)");
    hh_pulls[ii]->Fit("gausFitt2", "Q", "", -5, 5);
    
    // What it should be 
    TF1 *gausFit = new TF1("gausFitt", gausFitFunc, -5, 5, 3);
    gausFit->SetParameters( gausFit2->GetParameter(0), 0, 1);
    gausFit->SetLineColorAlpha( kRed, 0.35);
    
    gausFit->Draw("same");
    gausFit2->Draw("same");
  }
  printPNGandPDF( c1,  pngName, imgoutdir );

  fileout->Write();
  fileout->Close();
  
 return; 

}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Bunch of junk methods that I'm still playing around with

//________________________________________________________________________________________________________________________________

Double_t gausFitFunc(Double_t *v, Double_t *par) {

  Double_t arg = (v[0] - par[1])/par[2];

  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
  return fitval;

}

//________________________________________________________________________________________________________________________________

Double_t calcSig( Double_t confCutL, Double_t plateau ){

  Double_t sig   = (1-confCutL) * plateau  ; 

  return sig;
}
//________________________________________________________________________________________________________________________________


Double_t calcTot( TH1D* confLevs, Double_t confCutL, Double_t plateau ){

  Int_t bmin = confLevs->FindBin( confCutL );
  Int_t bmax = confLevs->FindBin( 1        );
  Double_t tot = confLevs->Integral( bmin, bmax, "width" );

  return tot;
}

//________________________________________________________________________________________________________________________________

Double_t calcNoise( TH1D* confLevs, Double_t confCutL, Double_t plateau ){
  
  Double_t sig = calcSig(confCutL, plateau);
  Double_t tot = calcTot(confLevs, confCutL, plateau);
  
  Double_t noise = tot - sig;

  return noise;
}


//________________________________________________________________________________________________________________________________

Double_t calcSNR( TH1D* confLevs, Double_t confCutL, Double_t plateau ){

  Double_t sig   = calcSig(confCutL, plateau);
  Double_t noise = calcNoise(confLevs, confCutL, plateau);


  Double_t snr = sig/noise;
  return snr;
}

//________________________________________________________________________________________________________________________________

Double_t calcSigPct( TH1D* confLevs, Double_t confCutL, Double_t plateau ){

  Double_t sig = calcSig(confCutL, plateau);
  Double_t tot = calcTot(confLevs, confCutL, plateau);

  Double_t sigpct = sig/tot * 100;
  return sigpct;
}


//______________________________________________________________________________
