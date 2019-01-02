#ifndef particlehistofiller_hh
#define particlehistofiller_hh

#include "ParticleHistoPlotter.h"
#include "KaonsHistoPlotter.h"
#include "TLorentzVector.h"

class ParticleHistoFiller{

 public:
  ParticleHistoFiller();
  ParticleHistoFiller(const char*);
  virtual ~ParticleHistoFiller();

  const char* temp_hname;

  ParticleHistoPlotter h_el_plotter;
  ParticleHistoPlotter h_pr_plotter;
  ParticleHistoPlotter h_kp_plotter;
  ParticleHistoPlotter h_km_plotter;
  ParticleHistoPlotter h_km0_plotter;
  ParticleHistoPlotter h_km3_plotter;
  KaonsHistoPlotter h_kaons_plotter;
  KaonsHistoPlotter h_kaons0_plotter;
  KaonsHistoPlotter h_kaons3_plotter;
  
  
 public:

  void initalizeAllHistos();

  void FillElectronHisto(TLorentzVector);
  void FillProtonHisto(TLorentzVector);
  void FillKaonPlusHisto(TLorentzVector);
  void FillKaonMinusHisto(TLorentzVector);
  void FillKaonMinus0Histo(TLorentzVector);
  void FillKaonMinus3Histo(TLorentzVector);
  void FillKaonPairHisto(std::vector<TLorentzVector>);
  void FillKaonPair0Histo(std::vector<TLorentzVector>);
  void FillKaonPair3Histo(std::vector<TLorentzVector>);
  
  void PrintAllHisto();
  
};
#endif
