#ifndef Analysis_h
#define Analysis_h

#include <cmath>
#include <iostream>
#include <string>
#include <TSpectrum.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TTree.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <algorithm>


// Header file for the classes stored in the TTree if any.

class Analysis {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           BB10Mul;
   Int_t           BB10Det[128];   //[BB10Mul]
   Int_t           BB10Strip[128];   //[BB10Mul]
   Int_t           BB10Channel[128];   //[BB10Mul]
   Int_t           BB10ADC[128];   //[BB10Mul]
   Float_t         BB10Energy[128];   //[BB10Mul]
   Int_t           QQQ5Mul;
   Char_t          QQQ5Upstream[128];   //[QQQ5Mul]
   Int_t           QQQ5Det[128];   //[QQQ5Mul]
   Int_t           QQQ5Ring[128];   //[QQQ5Mul]
   Int_t           QQQ5RingChannel[128];   //[QQQ5Mul]
   Int_t           QQQ5Sector[128];   //[QQQ5Mul]
   Int_t           QQQ5SectorChannel[128];   //[QQQ5Mul]
   Int_t           QQQ5RingADC[128];   //[QQQ5Mul]
   Float_t         QQQ5RingEnergy[128];   //[QQQ5Mul]
   Int_t           QQQ5SectorADC[128];   //[QQQ5Mul]
   Float_t         QQQ5SectorEnergy[128];   //[QQQ5Mul]
   Float_t         QQQ5Angle[128];   //[QQQ5Mul]
   Int_t           SX3Mul;
   Char_t          SX3Upstream[128];   //[SX3Mul]
   Int_t           SX3Det[128];   //[SX3Mul]
   Int_t           SX3Sector[128];   //[SX3Mul]
   Int_t           SX3SectorChannel[128];   //[SX3Mul]
   Int_t           SX3SectorADC[128];   //[SX3Mul]
   Float_t         SX3SectorEnergy[128];   //[SX3Mul]
   Int_t           SX3Strip[128];   //[SX3Mul]
   Int_t           SX3StripLeftChannel[128];   //[SX3Mul]
   Int_t           SX3StripRightChannel[128];   //[SX3Mul]
   Int_t           SX3StripLeftADC[128];   //[SX3Mul]
   Int_t           SX3StripRightADC[128];   //[SX3Mul]
   Int_t         SX3StripEnergy[128];   //[SX3Mul]
   Int_t           icdE;
   Int_t           icE;
   Int_t           icWireX;
   Int_t           icWireY;
   Float_t         icPositionX;
   Float_t         icPositionY;
   Float_t         icPositionWeightedX;
   Float_t         icPositionWeightedY;
   Int_t           tdcIC;
   Int_t           tdcGRETINA;
   Int_t           tdcRF;
   Int_t           tdcSilicon;
   ULong64_t       timeStamp;

   // List of branches
   TBranch        *b_BB10Mul;   //!
   TBranch        *b_BB10Det;   //!
   TBranch        *b_BB10Strip;   //!
   TBranch        *b_BB10Channel;   //!
   TBranch        *b_BB10ADC;   //!
   TBranch        *b_BB10Energy;   //!
   TBranch        *b_QQQ5Mul;   //!
   TBranch        *b_QQQ5Upstream;   //!
   TBranch        *b_QQQ5Det;   //!
   TBranch        *b_QQQ5Ring;   //!
   TBranch        *b_QQQ5RingChannel;   //!
   TBranch        *b_QQQ5Sector;   //!
   TBranch        *b_QQQ5SectorChannel;   //!
   TBranch        *b_QQQ5RingADC;   //!
   TBranch        *b_QQQ5RingEnergy;   //!
   TBranch        *b_QQQ5SectorADC;   //!
   TBranch        *b_QQQ5SectorEnergy;   //!
   TBranch        *b_QQQ5Angle;   //!
   TBranch        *b_SX3Mul;   //!
   TBranch        *b_SX3Upstream;   //!
   TBranch        *b_SX3Det;   //!
   TBranch        *b_SX3Sector;   //!
   TBranch        *b_SX3SectorChannel;   //!
   TBranch        *b_SX3SectorADC;   //!
   TBranch        *b_SX3SectorEnergy;   //!
   TBranch        *b_SX3Strip;   //!
   TBranch        *b_SX3StripLeftChannel;   //!
   TBranch        *b_SX3StripRightChannel;   //!
   TBranch        *b_SX3StripLeftADC;   //!
   TBranch        *b_SX3StripRightADC;   //!
   TBranch        *b_SX3StripEnergy;   //!
   TBranch        *b_icdE;   //!
   TBranch        *b_icE;   //!
   TBranch        *b_icWireX;   //!
   TBranch        *b_icWireY;   //!
   TBranch        *b_icPositionX;   //!
   TBranch        *b_icPositionY;   //!
   TBranch        *b_icPositionWeightedX;   //!
   TBranch        *b_icPositionWeightedY;   //!
   TBranch        *b_tdcIC;   //!
   TBranch        *b_tdcGRETINA;   //!
   TBranch        *b_tdcRF;   //!
   TBranch        *b_tdcSilicon;   //!
   TBranch        *b_timeStamp;   //!

   Analysis(TTree *tree=0);
   virtual ~Analysis();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);


// Variables
private:
    // Cuts
    TCutG* icdEECut;
};

#endif

#ifdef Analysis_cxx
Analysis::Analysis(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Run0073.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Run0073.root");
      }
      f->GetObject("data",tree);

   }
   Init(tree);
}

Analysis::~Analysis()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Analysis::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Analysis::LoadTree(Long64_t entry)
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

void Analysis::Init(TTree *tree)
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

   fChain->SetBranchAddress("BB10Mul", &BB10Mul, &b_BB10Mul);
   fChain->SetBranchAddress("BB10Det", BB10Det, &b_BB10Det);
   fChain->SetBranchAddress("BB10Strip", BB10Strip, &b_BB10Strip);
   fChain->SetBranchAddress("BB10Channel", BB10Channel, &b_BB10Channel);
   fChain->SetBranchAddress("BB10ADC", BB10ADC, &b_BB10ADC);
   fChain->SetBranchAddress("BB10Energy", BB10Energy, &b_BB10Energy);
   fChain->SetBranchAddress("QQQ5Mul", &QQQ5Mul, &b_QQQ5Mul);
   fChain->SetBranchAddress("QQQ5Upstream", QQQ5Upstream, &b_QQQ5Upstream);
   fChain->SetBranchAddress("QQQ5Det", QQQ5Det, &b_QQQ5Det);
   fChain->SetBranchAddress("QQQ5Ring", QQQ5Ring, &b_QQQ5Ring);
   fChain->SetBranchAddress("QQQ5RingChannel", QQQ5RingChannel, &b_QQQ5RingChannel);
   fChain->SetBranchAddress("QQQ5Sector", QQQ5Sector, &b_QQQ5Sector);
   fChain->SetBranchAddress("QQQ5SectorChannel", QQQ5SectorChannel, &b_QQQ5SectorChannel);
   fChain->SetBranchAddress("QQQ5RingADC", QQQ5RingADC, &b_QQQ5RingADC);
   fChain->SetBranchAddress("QQQ5RingEnergy", QQQ5RingEnergy, &b_QQQ5RingEnergy);
   fChain->SetBranchAddress("QQQ5SectorADC", QQQ5SectorADC, &b_QQQ5SectorADC);
   fChain->SetBranchAddress("QQQ5SectorEnergy", QQQ5SectorEnergy, &b_QQQ5SectorEnergy);
   fChain->SetBranchAddress("QQQ5Angle", QQQ5Angle, &b_QQQ5Angle);
   fChain->SetBranchAddress("SX3Mul", &SX3Mul, &b_SX3Mul);
   fChain->SetBranchAddress("SX3Upstream", SX3Upstream, &b_SX3Upstream);
   fChain->SetBranchAddress("SX3Det", SX3Det, &b_SX3Det);
   fChain->SetBranchAddress("SX3Sector", SX3Sector, &b_SX3Sector);
   fChain->SetBranchAddress("SX3SectorChannel", SX3SectorChannel, &b_SX3SectorChannel);
   fChain->SetBranchAddress("SX3SectorADC", SX3SectorADC, &b_SX3SectorADC);
   fChain->SetBranchAddress("SX3SectorEnergy", SX3SectorEnergy, &b_SX3SectorEnergy);
   fChain->SetBranchAddress("SX3Strip", SX3Strip, &b_SX3Strip);
   fChain->SetBranchAddress("SX3StripLeftChannel", SX3StripLeftChannel, &b_SX3StripLeftChannel);
   fChain->SetBranchAddress("SX3StripRightChannel", SX3StripRightChannel, &b_SX3StripRightChannel);
   fChain->SetBranchAddress("SX3StripLeftADC", SX3StripLeftADC, &b_SX3StripLeftADC);
   fChain->SetBranchAddress("SX3StripRightADC", SX3StripRightADC, &b_SX3StripRightADC);
   fChain->SetBranchAddress("SX3StripEnergy", SX3StripEnergy, &b_SX3StripEnergy);
   fChain->SetBranchAddress("icdE", &icdE, &b_icdE);
   fChain->SetBranchAddress("icE", &icE, &b_icE);
   fChain->SetBranchAddress("icWireX", &icWireX, &b_icWireX);
   fChain->SetBranchAddress("icWireY", &icWireY, &b_icWireY);
   fChain->SetBranchAddress("icPositionX", &icPositionX, &b_icPositionX);
   fChain->SetBranchAddress("icPositionY", &icPositionY, &b_icPositionY);
   fChain->SetBranchAddress("icPositionWeightedX", &icPositionWeightedX, &b_icPositionWeightedX);
   fChain->SetBranchAddress("icPositionWeightedY", &icPositionWeightedY, &b_icPositionWeightedY);
   fChain->SetBranchAddress("tdcIC", &tdcIC, &b_tdcIC);
   fChain->SetBranchAddress("tdcGRETINA", &tdcGRETINA, &b_tdcGRETINA);
   fChain->SetBranchAddress("tdcRF", &tdcRF, &b_tdcRF);
   fChain->SetBranchAddress("tdcSilicon", &tdcSilicon, &b_tdcSilicon);
   fChain->SetBranchAddress("timeStamp", &timeStamp, &b_timeStamp);
   Notify();
}

Bool_t Analysis::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Analysis::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Analysis::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Analysis_cxx
