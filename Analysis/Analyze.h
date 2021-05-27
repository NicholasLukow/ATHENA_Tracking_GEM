//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu May 20 14:11:13 2021 by ROOT version 6.22/02
// from TTree tracks/FastSim Eval => tracks
// found on file: out_TrackingStudy_Beast_FastSimEval.root
//////////////////////////////////////////////////////////

#ifndef Analyze_h
#define Analyze_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Analyze {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           gtrackID;
   Int_t           gflavor;
   Float_t         gpx;
   Float_t         gpy;
   Float_t         gpz;
   Float_t         gvx;
   Float_t         gvy;
   Float_t         gvz;
   Float_t         gvt;
   Int_t           trackID;
   Int_t           charge;
   Int_t           nhits;
   Float_t         px;
   Float_t         py;
   Float_t         pz;
   Float_t         pcax;
   Float_t         pcay;
   Float_t         pcaz;
   Float_t         dca2d;
   Int_t           nHit_G4HIT_BMT;
   Int_t           nHit_G4HIT_SVTX;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_gtrackID;   //!
   TBranch        *b_gflavor;   //!
   TBranch        *b_gpx;   //!
   TBranch        *b_gpy;   //!
   TBranch        *b_gpz;   //!
   TBranch        *b_gvx;   //!
   TBranch        *b_gvy;   //!
   TBranch        *b_gvz;   //!
   TBranch        *b_gvt;   //!
   TBranch        *b_trackID;   //!
   TBranch        *b_charge;   //!
   TBranch        *b_nhits;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_pcax;   //!
   TBranch        *b_pcay;   //!
   TBranch        *b_pcaz;   //!
   TBranch        *b_dca2d;   //!
   TBranch        *b_nHit_G4HIT_BMT;   //!
   TBranch        *b_nHit_G4HIT_SVTX;   //!

   Analyze(TTree *tree=0);
   virtual ~Analyze();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Analyze_cxx
Analyze::Analyze(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("out_TrackingStudy_Beast_FastSimEval.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("out_TrackingStudy_Beast_FastSimEval.root");
      }
      f->GetObject("tracks",tree);

   }
   Init(tree);
}

Analyze::~Analyze()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analyze::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analyze::LoadTree(Long64_t entry)
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

void Analyze::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("gtrackID", &gtrackID, &b_gtrackID);
   fChain->SetBranchAddress("gflavor", &gflavor, &b_gflavor);
   fChain->SetBranchAddress("gpx", &gpx, &b_gpx);
   fChain->SetBranchAddress("gpy", &gpy, &b_gpy);
   fChain->SetBranchAddress("gpz", &gpz, &b_gpz);
   fChain->SetBranchAddress("gvx", &gvx, &b_gvx);
   fChain->SetBranchAddress("gvy", &gvy, &b_gvy);
   fChain->SetBranchAddress("gvz", &gvz, &b_gvz);
   fChain->SetBranchAddress("gvt", &gvt, &b_gvt);
   fChain->SetBranchAddress("trackID", &trackID, &b_trackID);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("nhits", &nhits, &b_nhits);
   fChain->SetBranchAddress("px", &px, &b_px);
   fChain->SetBranchAddress("py", &py, &b_py);
   fChain->SetBranchAddress("pz", &pz, &b_pz);
   fChain->SetBranchAddress("pcax", &pcax, &b_pcax);
   fChain->SetBranchAddress("pcay", &pcay, &b_pcay);
   fChain->SetBranchAddress("pcaz", &pcaz, &b_pcaz);
   fChain->SetBranchAddress("dca2d", &dca2d, &b_dca2d);
   fChain->SetBranchAddress("nHit_G4HIT_BMT", &nHit_G4HIT_BMT, &b_nHit_G4HIT_BMT);
   fChain->SetBranchAddress("nHit_G4HIT_SVTX", &nHit_G4HIT_SVTX, &b_nHit_G4HIT_SVTX);
   Notify();
}

Bool_t Analyze::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analyze::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analyze::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analyze_cxx
