#include "PhysicsHistoPlotter.h"

#include "TH1D.h"
#include "TH2D.h"

PhysicsHistoPlotter::PhysicsHistoPlotter(){


}

PhysicsHistoPlotter::~PhysicsHistoPlotter(){

}


void PhysicsHistoPlotter::InitPhysicsHistos( const char* temp_title ){

  f = new TFile(Form("h_Xs_%s.root",temp_title),"RECREATE");
  
  h_eh1X = new TH1D(Form("h_eh1X_%s",temp_title),Form("h_eh1X_%s",temp_title), 200, 0.0, 0.0 );
  h_eh2X = new TH1D(Form("h_eh2X_%s",temp_title),Form("h_eh2X_%s",temp_title), 200, 0.0, 0.0 );
  
}


void PhysicsHistoPlotter::DisplayHistos(const char* temp_title ){

  TCanvas *c1 = new TCanvas(Form("c_Xs_%s",temp_title),Form("c_Xs_%s",temp_title), 900, 900 );
  h_eh1X->Draw();
  h_eh2X->Draw();

  f->Write();
  f->Save();
  f->Close();
   

}
