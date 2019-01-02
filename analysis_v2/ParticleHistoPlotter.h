#ifndef particlehistoplotter_hh
#define particlehistoplotter_hh

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

class ParticleHistoPlotter{

 public:
  //std::string part_name;
  ParticleHistoPlotter();
  virtual ~ParticleHistoPlotter();
  
  TH1D *h_p;
  TH1D *h_theta;
  TH1D *h_phi;

  TFile *f;
  
  TH2D *h2_ptheta;
  TH2D *h2_pphi;
  TH2D *h2_thetaphi;

 public:
  void initializeHistos(const char*);
  void InitializeKaonPairHistos(const char*);
  void RescaleAxis();
  void DisplayHistos(const char*);
  
};
#endif
