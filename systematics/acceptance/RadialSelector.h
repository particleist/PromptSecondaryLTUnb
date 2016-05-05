//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Apr 30 17:01:38 2016 by ROOT version 6.06/02
// from TTree DecayTree/DecayTree
// found on file: /data/Phys/data/mcd02kpi_tracks9_merged.root
//////////////////////////////////////////////////////////

#ifndef RadialSelector_h
#define RadialSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TH1D.h>
#include <TCanvas.h>

class RadialSelector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Int_t> D_ID = {fReader, "D_ID"};
   TTreeReaderValue<Double_t> D_PX = {fReader, "D_PX"};
   TTreeReaderValue<Double_t> D_PY = {fReader, "D_PY"};
   TTreeReaderValue<Double_t> D_PZ = {fReader, "D_PZ"};
   TTreeReaderValue<Double_t> D_PT = {fReader, "D_PT"};
   TTreeReaderValue<Double_t> D_E = {fReader, "D_E"};
   TTreeReaderValue<Double_t> D_ETA = {fReader, "D_ETA"};
   TTreeReaderValue<Double_t> D_PHI = {fReader, "D_PHI"};
   TTreeReaderValue<Double_t> D_MM = {fReader, "D_MM"};
   TTreeReaderValue<Double_t> D_VX = {fReader, "D_VX"};
   TTreeReaderValue<Double_t> D_VY = {fReader, "D_VY"};
   TTreeReaderValue<Double_t> D_VZ = {fReader, "D_VZ"};
   TTreeReaderValue<Double_t> D_PVX = {fReader, "D_PVX"};
   TTreeReaderValue<Double_t> D_PVY = {fReader, "D_PVY"};
   TTreeReaderValue<Double_t> D_PVZ = {fReader, "D_PVZ"};

   TTreeReaderValue<Int_t> D_FROMB = {fReader, "D_FROMB"};
   TTreeReaderValue<Double_t> D_BPVDIRA = {fReader, "D_BPVDIRA"};
   TTreeReaderValue<Double_t> D_BPVLTIME = {fReader, "D_BPVLTIME"};
   // Decay time w.r.t. BPV, t = (flight distance) * (mass) / (momentum).

   // Internal counter
   double_t mmTotal;
   int      eventCount;
   TCanvas *c1;
   TH1D *histMM;
   TH1D *histPVz;
   TH1D *histPVr;
   TH1D *histEVz;
   TH1D *histEVr;
   TH1D *histPEta;
   TH1D *histEVz_FromB;
   TH1D *histEVr_FromB;
   TH1D *histPEta_FromB;
   TH1D *histLifetime;
   TH1D *histAcceptance;
   TH1D *histAcceptanceV;
   TH1D *histAcceptanceRatio;
   TH1D *histCount;


   // Some constants
   const Double_t HIST_ACCEPTANCE_MAX = 0.5;
   
 RadialSelector(TTree * /*tree*/ =0) : mmTotal(0), eventCount(0),
     histMM(nullptr),
     histPVz(nullptr),
     histPVr(nullptr),
     histEVz(nullptr),
     histEVr(nullptr),
     histPEta(nullptr),
     histEVz_FromB(nullptr),
     histEVr_FromB(nullptr),
     histPEta_FromB(nullptr),
     histLifetime(nullptr),
     histAcceptance(nullptr),
     histAcceptanceV(nullptr),
     histAcceptanceRatio(nullptr),
     histCount(nullptr)
     { }
   virtual ~RadialSelector() {}
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   // Util methods
   Double_t GetRadius(Double_t z, Double_t px, Double_t py, Double_t pz, Double_t vx, Double_t vy, Double_t vz);
   Double_t GetZforRadius(Double_t R, Double_t px, Double_t py, Double_t pz, Double_t vx, Double_t vy, Double_t vz);
      
   ClassDef(RadialSelector,0);
};

#endif

#ifdef RadialSelector_cxx
void RadialSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t RadialSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef RadialSelector_cxx
