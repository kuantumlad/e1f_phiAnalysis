#ifndef physicshistoplotter_hh
#define physicshistoplotter_hh

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

class PhysicsHistoPlotter{


 public:
  PhysicsHistoPlotter();
  virtual ~PhysicsHistoPlotter();

  TH1D *h_eh1X;
  TH1D *h_eh2X;

  TFile *f;

  
 public:
  void InitPhysicsHistos( const char* );
  void DisplayHistos( const char* );
    
  
};
#endif
