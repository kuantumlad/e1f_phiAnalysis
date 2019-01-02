#ifndef ApplyKinFit_h
#define ApplyKinFit_h

#include <iostream>
#include <vector>

#include "ExclusiveEvents.C"
#include "CLAS_KinFitter/include/KinFitter.h"
#include "CLAS_KinFitter/include/res-errs_clas.h"
#include "TLorentzVector.h"

std::vector<TLorentzVector> ApplyKinFit(TLorentzVector lv_el, TLorentzVector lv_pr, TLorentzVector lv_kp, TLorentzVector lv_km, int topology, KinFitter *kinFit ){

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

  TLorentzVector lvkf_e = kinFit->KFParticlesOut[0].P;
  TLorentzVector lvkf_pr = kinFit->KFParticlesOut[1].P;
  TLorentzVector lvkf_kp = kinFit->KFParticlesOut[2].P;
  TLorentzVector lvkf_km = kinFit->KFParticlesOut[3].P;
  TLorentzVector lvkf_phi = lvkf_kp + lvkf_km;
  TLorentzVector lvkf_ekpX = lv_ebeam + lv_target - lvkf_e - lvkf_kp;

  std::vector<TLorentzVector> *lv_final;
  lv_final.push_back(lvkf_e);
  lv_final.push_back(lvkf_pr);
  lv_final.push_back(lvkf_kp);
  lv_final.push_back(lvkf_km);

      
}
#endif
