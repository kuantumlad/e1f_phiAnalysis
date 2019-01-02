#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TFile.h>
#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <TString.h>
#include "TROOT.h"
#include "TBrowser.h"
#include <iostream>
#include <string>
#include "TChain.h"
#include "TRegexp.h"
#include "TObject.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <map>
#include "TCanvas.h"
#include "TPad.h"
#include "TBox.h"
////////////////////////////
////USE FOR APPLICATION OF WEIGHTS
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
////////////////////////////

//////////////////////////
//USE TO TRAIN AND TEST MODEL
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
/////////////////////////

void e1fNNAnalysis( ){

  std::cout << " >> STARTING NN ANALYSIS " << std::endl;

  //LOAD SIGNAL AND BACKGROUND TREES

  TString outfileName("TMVA.root");
  TFile *outputFile = new TFile(outfileName,"RECREATE");
  
  TMVA::Factory *factory = new TMVA::Factory("TMVAClassification",outputFile,"!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G;D:AnalysisType=Classification");

  TMVA::DataLoader *dataloader = new TMVA::DataLoader("dataset");
  

}
