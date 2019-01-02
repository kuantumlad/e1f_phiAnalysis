#include "PhysicsHistoFiller.h"
#include "TH1D.h"
#include "TH2D.h"

PhysicsHistoFiller::PhysicsHistoFiller(){


}

PhysicsHistoFiller::PhysicsHistoFiller(const char* temp_name){

  temp_hname = temp_name;

}

PhysicsHistoFiller::~PhysicsHistoFiller(){


}


void PhysicsHistoFiller::InitAllHistos(){

  h_epX.InitPhysicsHistos(Form("epX_%s",temp_hname));
  h_epkpX.InitPhysicsHistos(Form("epkpX_%s", temp_hname));

}


void PhysicsHistoFiller::FillPhysics(std::vector<TLorentzVector> v_lv ){


}

void PhysicsHistoFiller::FillPhysicsH1( TLorentzVector temp_lv ){

  h_epX.h_eh1X->Fill( temp_lv.M() );

}

void PhysicsHistoFiller::FillPhysicsH2( TLorentzVector temp_lv ){

  h_epX.h_eh2X->Fill( temp_lv.M() );

}


void PhysicsHistoFiller::PrintAllHisto(){

  h_epX.DisplayHistos(Form("epX_%s",temp_hname));
  h_epkpX.DisplayHistos(Form("epkpX_%s",temp_hname));
  

}
