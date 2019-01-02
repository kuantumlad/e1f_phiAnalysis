//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 15 17:26:22 2018 by ROOT version 6.06/04
// from TTree events/events
// found on file: events.root
//////////////////////////////////////////////////////////

#ifndef ExclusiveEvents_h
#define ExclusiveEvents_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TLorentzVector.h"

class ExclusiveEvents {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           helicity;
   Int_t           topology;
   TLorentzVector  *electron;
   TLorentzVector  *proton;
   TLorentzVector  *kaon_pos;
   TLorentzVector  *kaon_neg;
   Float_t         alpha_kaon_pos;
   Float_t         alpha_kaon_neg;
   Float_t         alpha_proton;
   Float_t         dist_ecsf;
   Float_t         dist_ec_edep;
   Float_t         dist_vz;
   Float_t         dist_cc_theta;
   Float_t         dist_dcr1;
   Float_t         dist_dcr3;
   Float_t         dist_ecu;
   Float_t         dist_ecv;
   Float_t         dist_ecw;
   Float_t         dist_cc;

   // List of branches
   TBranch        *b_helicity;   //!
   TBranch        *b_topology;   //!
   TBranch        *b_electron;   //!
   TBranch        *b_proton;   //!
   TBranch        *b_kaon_pos;   //!
   TBranch        *b_kaon_neg;   //!
   TBranch        *b_alpha_kaon_pos;   //!
   TBranch        *b_alpha_kaon_neg;   //!
   TBranch        *b_alpha_proton;   //!
   TBranch        *b_dist_ecsf;   //!
   TBranch        *b_dist_ec_edep;   //!
   TBranch        *b_dist_vz;   //!
   TBranch        *b_dist_cc_theta;   //!
   TBranch        *b_dist_dcr1;   //!
   TBranch        *b_dist_dcr3;   //!
   TBranch        *b_dist_ecu;   //!
   TBranch        *b_dist_ecv;   //!
   TBranch        *b_dist_ecw;   //!
   TBranch        *b_dist_cc;   //!

   ExclusiveEvents(TTree *tree=0);
   virtual ~ExclusiveEvents();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ExclusiveEvents_cxx
ExclusiveEvents::ExclusiveEvents(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("events.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("events.root");
      }
      f->GetObject("events",tree);

   }
   Init(tree);
}

ExclusiveEvents::~ExclusiveEvents()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ExclusiveEvents::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ExclusiveEvents::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ExclusiveEvents::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   electron = 0;
   proton = 0;
   kaon_pos = 0;
   kaon_neg = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("helicity", &helicity, &b_helicity);
   fChain->SetBranchAddress("topology", &topology, &b_topology);
   fChain->SetBranchAddress("electron", &electron, &b_electron);
   fChain->SetBranchAddress("proton", &proton, &b_proton);
   fChain->SetBranchAddress("kaon_pos", &kaon_pos, &b_kaon_pos);
   fChain->SetBranchAddress("kaon_neg", &kaon_neg, &b_kaon_neg);
   fChain->SetBranchAddress("alpha_kaon_pos", &alpha_kaon_pos, &b_alpha_kaon_pos);
   fChain->SetBranchAddress("alpha_kaon_neg", &alpha_kaon_neg, &b_alpha_kaon_neg);
   fChain->SetBranchAddress("alpha_proton", &alpha_proton, &b_alpha_proton);
   fChain->SetBranchAddress("dist_ecsf", &dist_ecsf, &b_dist_ecsf);
   fChain->SetBranchAddress("dist_ec_edep", &dist_ec_edep, &b_dist_ec_edep);
   fChain->SetBranchAddress("dist_vz", &dist_vz, &b_dist_vz);
   fChain->SetBranchAddress("dist_cc_theta", &dist_cc_theta, &b_dist_cc_theta);
   fChain->SetBranchAddress("dist_dcr1", &dist_dcr1, &b_dist_dcr1);
   fChain->SetBranchAddress("dist_dcr3", &dist_dcr3, &b_dist_dcr3);
   fChain->SetBranchAddress("dist_ecu", &dist_ecu, &b_dist_ecu);
   fChain->SetBranchAddress("dist_ecv", &dist_ecv, &b_dist_ecv);
   fChain->SetBranchAddress("dist_ecw", &dist_ecw, &b_dist_ecw);
   fChain->SetBranchAddress("dist_cc", &dist_cc, &b_dist_cc);
   Notify();
}

Bool_t ExclusiveEvents::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ExclusiveEvents::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ExclusiveEvents::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ExclusiveEvents_cxx
