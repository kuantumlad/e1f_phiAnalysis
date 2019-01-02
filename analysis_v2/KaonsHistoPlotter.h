#ifndef kaonshistoplotter_hh
#define kaonshistoplotter_hh

#include "TFile.h"

class KaonsHistoPlotter{


 public:
  KaonsHistoPlotter();
  virtual ~KaonsHistoPlotter();

  TFile *f;
  TH1D *h_decayang;
  TH2D *h_decayangp;
  TH2D *h_decayangp_prkm;
  TH2D *h_decayang_prkm_kpkm;

 public:
  void InitKaonHistos(const char* );
  void DisplayHistos( const char* );
  
};
#endif
