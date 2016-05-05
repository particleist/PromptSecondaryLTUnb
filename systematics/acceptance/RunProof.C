#include <iostream>
#include <string>

#include "TFile.h"
#include "TDSet.h"
#include "TChain.h"
#include "TProof.h"

void RunProof() {
  
  std::string filename = "/data/Phys/data/mcd02kpi_tracks9_merged.root";
  bool fileNotFound = gSystem->AccessPathName(filename.c_str());
  // if local file is not found use the EOS version
  if (fileNotFound) {
    filename = "root://eoslhcb.cern.ch//eos/lhcb/user/m/malexand/d2hh/secondariesTagging/mcd02kpi_tracks9_merged.root";
  }
  
  TChain chain("TrackFilter/DecayTree");
  chain.Add(filename.c_str());
  TDSet dset(chain);

  auto p = TProof::Open("");
  p->Process(&dset, "RadialSelector.C");

}
