#define RadialSelector_cxx

// The class definition in RadialSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("RadialSelector.C")
// root> T->Process("RadialSelector.C","some options")
// root> T->Process("RadialSelector.C+")
//


#include "RadialSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <cmath>

void RadialSelector::Begin(TTree * /*tree*/)
{
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
}

void RadialSelector::SlaveBegin(TTree * /*tree*/)
{
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).

  TString option = GetOption();
  
  // The mass histogram
  histMM = new TH1D("DMassHist", "D Mass Histogram",
		    200, 1800, 1920);
  histMM->Sumw2();
  GetOutputList()->Add(histMM);

  // Histogram for the PV position
  histPVz = new TH1D("PVzHist", "PV z Histogram",
		     1000, -300, 300);
  histPVr = new TH1D("PVrHist", "PV r Histogram",
		     200, 0.0, 1.0);
  GetOutputList()->Add(histPVz);
  GetOutputList()->Add(histPVr);


  // Histogram for the end vertex position
  histEVz = new TH1D("EVzHist", "EV z Histogram Prompt D",
		     1000, -200, 200);
  histEVr = new TH1D("EVrHist", "EV r Histogram Prompt D",
		     1000, 0.0, 2.0);
  GetOutputList()->Add(histEVz);
  GetOutputList()->Add(histEVr);

  // histogram for the momentum distribution in ETA
  histPEta = new TH1D("PEtaHist", "Momentum Eta distribution Prompt D",
		     400, -0.1, 30);
  GetOutputList()->Add(histPEta);

  // Histogram for the end vertex position
  histEVz_FromB = new TH1D("EVzHist_FromB", "EV z Histogram from B",
		     1000, -200, 200);
  histEVr_FromB = new TH1D("EVrHist_FromB", "EV r Histogram from B",
		     1000, 0.0, 2.0);
  GetOutputList()->Add(histEVz_FromB);
  GetOutputList()->Add(histEVr_FromB);

  // histogram for the momentum distribution in ETA
  histPEta_FromB = new TH1D("PEtaHist_FromB", "Momentum Eta distribution from B",
		     400, -0.1, 30);
  GetOutputList()->Add(histPEta_FromB);

  // histogram for the lifetime
  histLifetime = new TH1D("LtimeHist", "Calculated lifetime / D_BPVLTIME",
		     500, 0, 2);
  GetOutputList()->Add(histLifetime);

  // histogram for the acceptance
  histAcceptance = new TH1D("AcceptanceHist", "Lifetime acceptance (extended following momentum direction)",
		     500, 0, HIST_ACCEPTANCE_MAX);
  GetOutputList()->Add(histAcceptance);

  // histogram for the acceptance ratio
  histAcceptanceRatio = new TH1D("AcceptanceRatioHist", "Ratio between max lifetime (following momentum) and actual value",
		     500, 0, 200);
  GetOutputList()->Add(histAcceptanceRatio);

  histAcceptanceV = new TH1D("AcceptanceVerticesHist", "Lifetime acceptance (extended following vertices direction)",
		     500, 0, HIST_ACCEPTANCE_MAX);
  GetOutputList()->Add(histAcceptanceV);

  // Counters
  histCount = new TH1D("histCount", "Counters", 1, 0, 1);
  GetOutputList()->Add(histCount);
}


/**
 * @see RadialSelector.h
 */
Double_t RadialSelector::GetZforRadius(Double_t R,
				       const TVector3 &pv,
				       const TVector3 &ev) {

  Double_t a = ev.Perp2();
  Double_t b = 2.0 * (pv.X() *ev.X() + pv.Y() * ev.Y());
  Double_t c = pv.Perp2() - R * R;
  Double_t delta = b * b - 4 * a * c;

  if (delta < 0) {
    return nan("");
  }

  // Only the larger root is needed, looking in the forward direction...
  Double_t Z1 = (- b + sqrt(delta)) / (2 * a);
  Double_t Z2 = (- b - sqrt(delta)) / (2 * a);
  return ev.Z() * max(Z1, Z2) + pv.Z();
}

/**
 * @see RadialSelector.h
 */
Double_t RadialSelector::GetRadius(Double_t z0,
				   const TVector3 &pv,
				   const TVector3 &ev) {
  Double_t tmp = (z0 - pv.Z()) / ev.Z();
  Double_t rad = sqrt(pow(pv.X() + tmp * ev.X(), 2)
		      + pow( pv.Y() + tmp * ev.Y(), 2));
  return rad;
}

Bool_t RadialSelector::Process(Long64_t entry)
{
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  
  // General stuff
  fReader.SetEntry(entry);

  // Using the bin at value 0 tocount the number of events
  // Ugly to use histograms for that but this is merged automatically
  // by Proof
  histCount->Fill(0);
  
  // Mass stuff
  histMM->Fill(*D_MM);

  // Constructing the vectors
  TVector3 p(*D_PX, *D_PY, *D_PZ);
  TVector3 primaryVertex(*D_PVX, *D_PVY, *D_PVZ);
  TVector3 endVertex(*D_VX, *D_VY, *D_VZ);
  TVector3 diffVertex = endVertex - primaryVertex;
  
  // PV stats
  histPVr->Fill(primaryVertex.Perp());
  histPVz->Fill(primaryVertex.Z());
  
  if (*D_FROMB) {
    // Now checking the end vertex position
    histEVr_FromB->Fill(endVertex.Perp());
    histEVz_FromB->Fill(endVertex.Z());
    histPEta_FromB->Fill(p.Eta());
  } else {
    // Now checking the end vertex position
    histEVr->Fill(endVertex.Perp());
    histEVz->Fill(endVertex.Z());
    histPEta->Fill(p.Eta());
  }

  // Checking the lifetime
  Double_t fd = diffVertex.Mag();
  const Double_t c = TMath::C() * 1e-6;// * 1e3 / 1e9;//We need mm/ns in LHCb Units
  Double_t ltime = fd * (*D_MM) / (p.Mag() * c);
  histLifetime->Fill(ltime / *D_BPVLTIME);

  // Now establishing the acceptance
  const Double_t radiusCut = 4.0; // Unit is mm
  auto  zAtCutFollowingP = GetZforRadius(radiusCut, primaryVertex, p);

  Double_t acceptanceRatio = (zAtCutFollowingP - primaryVertex.Z()) / (endVertex.Z() - primaryVertex.Z());
  Double_t acceptance = *D_BPVLTIME * acceptanceRatio;
   
  histAcceptance->Fill(acceptance);
  histAcceptanceRatio->Fill(acceptanceRatio);
  
  // Checking the difference if we use the vertices direction instead
  auto  zAtCutFollowingVertices = GetZforRadius(radiusCut, primaryVertex, diffVertex);
  histAcceptanceV->Fill(*D_BPVLTIME * (zAtCutFollowingVertices - primaryVertex.Z()) / (endVertex.Z() - primaryVertex.Z()));

  // And we're done...
  return kTRUE;
}

void RadialSelector::SlaveTerminate()
{
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.
  
}

void RadialSelector::Terminate()
{
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.

  // Printing the total number of events (bin 1).
  std::cout << "Hist event count: " << histCount->GetBinContent(1) << std::endl;
  
  TCanvas *c1 = new TCanvas("c1","Histogram checks",200,10,700,900);
  c1->Divide(2,3);
  c1->cd(1);
  histMM->Fit("gaus");
  histMM->Draw();
  c1->cd(3);
  histPVz->Fit("gaus");
  histPVz->Draw();
  c1->cd(4);
  histPVr->Draw();
  c1->cd(5);
  histEVz->SetStats(0);
  histEVz->DrawNormalized();
  histEVz_FromB->SetStats(0);
  histEVz_FromB->SetLineColor(kGreen);
  histEVz_FromB->DrawNormalized("SAME");
  auto legendEVz = new TLegend(0.5,0.8,0.89,0.89);
  legendEVz->SetFillColor(0);
  legendEVz->AddEntry(histEVz);
  legendEVz->AddEntry(histEVz_FromB);
  legendEVz->Draw("SAME");

  c1->cd(6);
  histEVr->SetStats(0);
  histEVr->DrawNormalized();
  histEVr_FromB->SetStats(0);
  histEVr_FromB->SetLineColor(kGreen);
  histEVr_FromB->DrawNormalized("SAME");
  auto legendEVr = new TLegend(0.5,0.8,0.89,0.89);
  legendEVr->SetFillColor(0);
  legendEVr->AddEntry(histEVr);
  legendEVr->AddEntry(histEVr_FromB);
  legendEVr->Draw("SAME");

  c1->cd(2);
  histPEta->SetStats(0);
  histPEta->DrawNormalized();
  histPEta_FromB->SetStats(0);
  histPEta_FromB->SetLineColor(kGreen);
  histPEta_FromB->DrawNormalized("SAME");
  auto legendPEta = new TLegend(0.5,0.8,0.89,0.89);
  legendPEta->SetFillColor(0);
  legendPEta->AddEntry(histPEta);
  legendPEta->AddEntry(histPEta_FromB);
  legendPEta->Draw("SAME");
  c1->Update();
  
  TCanvas *c2 = new TCanvas("c2","Lifetime check",200,10,700,900);
  c2->SetLogy();
  histLifetime->DrawNormalized();
  c2->Update();
  
  TCanvas *c3 = new TCanvas("c3","Acceptance",200,10,700,900);
  c3->SetLogy();
  c3->SetTitle("Lifetime acceptance (ns)");
  histAcceptance->SetStats(0);
  histAcceptance->DrawNormalized();
  histAcceptanceV->SetStats(0);
  histAcceptanceV->SetLineColor(kGreen);
  histAcceptanceV->DrawNormalized("SAME");
  auto legendLf = new TLegend(0.3,0.85,0.89,0.89);
  legendLf->SetFillColor(0);
  legendLf->AddEntry(histAcceptance);
  legendLf->AddEntry(histAcceptanceV);
  legendLf->Draw("SAME");
  c3->Update();

  TCanvas *c4 = new TCanvas("c4","Ratio between lifetime and max lifetime",200,10,700,900);
  c4->SetLogy();
  histAcceptanceRatio->DrawNormalized();
  c4->Update();

}
