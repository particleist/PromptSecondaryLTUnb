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

TH1D* RadialSelector::NewHist(const char *key,
			      const char *description,
			      Int_t nbbins,
			      Double_t min,
			      Double_t max) {
  TH1D *h = new TH1D(key, description, nbbins, min, max);
  GetOutputList()->Add(h);
  return h;
}

TH1D* RadialSelector::H(const char *key) {
  // Getting our histogram from the output list
  TList *olist = GetOutputList();
  TIter next(olist);
  while (TObject *obj = next()) {
    if (strcmp(obj->GetName(), key) == 0) {
      return (TH1D *)obj;
    }
  }
  return nullptr;
}

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
  auto histMM = NewHist("MM", "D Mass Histogram", 200, 1800, 1920);
  histMM->Sumw2();

  // Histogram for the PV position
  NewHist("PVz", "PV z Histogram", 1000, -300, 300);
  NewHist("PVr", "PV r Histogram", 200, 0.0, 1.0);

  // Histogram for the end vertex position
  NewHist("EVz", "EV z Histogram Prompt D", 1000, -200, 200);
  NewHist("EVr", "EV r Histogram Prompt D", 1000, 0.0, 2.0);

  // histogram for the momentum distribution in ETA
  NewHist("PEta", "Momentum Eta distribution Prompt D", 400, -0.1, 30);

  // Histogram for the end vertex position
  NewHist("EVz_FromB", "EV z Histogram from B", 1000, -200, 200);
  NewHist("EVr_FromB", "EV r Histogram from B", 1000, 0.0, 2.0);

  // histogram for the momentum distribution in ETA
  NewHist("PEta_FromB", "Momentum Eta distribution from B", 400, -0.1, 30);

  // histogram for the lifetime
  NewHist("Lifetime", "Calculated lifetime / D_BPVLTIME", 500, 0, 2);

  // histogram for the acceptance
  NewHist("Acceptance", "Lifetime acceptance (extended following momentum direction)",
	  500, 0, HIST_ACCEPTANCE_MAX);

  // histogram for the acceptance lower lifetime
  NewHist("AcceptanceZ", "Lifetime acceptance (extended following momentum direction)",
	  500, 0, HIST_ACCEPTANCE_MAXZOOM);
  
  // histogram for the acceptance ratio
  NewHist("AcceptanceRatio", "Ratio between max lifetime (following momentum) and actual value",
	  500, 0, 200);
  
  NewHist("AcceptanceV", "Lifetime acceptance (extended following vertices direction)",
	  500, 0, HIST_ACCEPTANCE_MAX);
  
  NewHist("AcceptanceVZ", "Lifetime acceptance (extended following vertices direction)",
	  500, 0, HIST_ACCEPTANCE_MAXZOOM);

  // Counters
  auto histCount = NewHist("Count", "Counters", 1, 0, 1);
  histCount->SetCanExtend(TH1::kAllAxes); 
}


/**
 * @see RadialSelector.h
 */
Double_t RadialSelector::GetZforRadius(Double_t R,
				       const TVector3 &point,
				       const TVector3 &direction) {

  Double_t a = direction.Perp2();
  Double_t b = 2.0 * (point.X() *direction.X() + point.Y() * direction.Y());
  Double_t c = point.Perp2() - R * R;
  Double_t delta = b * b - 4 * a * c;

  if (delta < 0) {
    return nan("");
  }

  Double_t Z1 = (- b + sqrt(delta)) / (2 * a);
  Double_t Z2 = (- b - sqrt(delta)) / (2 * a);
  Double_t Zsel1, Zsel2, Z;
  Zsel1 = direction.Z() * Z1 + point.Z();
  Zsel2 = direction.Z() * Z2 + point.Z();
  if (direction.Z() < 0) {
    Z = min(Zsel1, Zsel2);
  } else {
    Z = max(Zsel1, Zsel2);
  }
  
  return Z;
}

/**
 * @see RadialSelector.h
 */
Double_t RadialSelector::GetRadius(Double_t z0,
				   const TVector3 &point,
				   const TVector3 &direction) {
  Double_t tmp = (z0 - point.Z()) / direction.Z();
  Double_t rad = sqrt(pow(point.X() + tmp * direction.X(), 2)
		      + pow( point.Y() + tmp * direction.Y(), 2));
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
  H("Count")->Fill("TREE_TOTAL", 1);  

  // Constructing the vectors
  TVector3 p(*D_PX, *D_PY, *D_PZ);
  TVector3 primaryVertex(*D_PVX, *D_PVY, *D_PVZ);
  TVector3 endVertex(*D_VX, *D_VY, *D_VZ);
  TVector3 diffVertex = endVertex - primaryVertex;
  
  // Now establishing the acceptance
  const Double_t radiusCut = 4.0; // Unit is mm
  bool getmin = ((endVertex.Z() - primaryVertex.Z()) < 0);
  auto  zAtCutFollowingP = GetZforRadius(radiusCut, primaryVertex, p);
  Double_t acceptanceRatio = (zAtCutFollowingP - primaryVertex.Z()) / (endVertex.Z() - primaryVertex.Z());
  Double_t acceptance = *D_BPVLTIME * acceptanceRatio;

  if (*D_BPVLTIME < 0) {
    H("Count")->Fill("NEG_TIME", 1);  
  }
  
  if (acceptanceRatio < 0) {
    H("Count")->Fill("NEG_ACCRATIO", 1);  
  }
  
  if (acceptanceRatio < 0 && *D_BPVLTIME <0) {
    H("Count")->Fill("NEG_BOTH", 1);  
  }

  if ((endVertex.Z() - primaryVertex.Z()) < 0) {
    H("Count")->Fill("EV_BEF_PV", 1);  
  }

  if ((endVertex.Z() - primaryVertex.Z()) < 0 && *D_BPVLTIME < 0) {
    H("Count")->Fill("NEG_TIME_EV_BEF_PV", 1);  
  }

  bool keepEvent= (*D_BPVLTIME > 0) && ((endVertex.Z() - primaryVertex.Z()) > 0);
  if (!keepEvent) {
    return kTRUE;
  }
    
  // Using the bin at value 0 tocount the number of events
  // Ugly to use histograms for that but this is merged automatically
  // by Proof
  H("Count")->Fill("TOTAL", 1);
  
  // Mass stuff
  H("MM")->Fill(*D_MM);
  
  // PV stats
  H("PVr")->Fill(primaryVertex.Perp());
  H("PVz")->Fill(primaryVertex.Z());
  
  if (*D_FROMB) {
    // Now checking the end vertex position
    H("EVr_FromB")->Fill(endVertex.Perp());
    H("EVz_FromB")->Fill(endVertex.Z());
    H("PEta_FromB")->Fill(p.Eta());
  } else {
    // Now checking the end vertex position
    H("EVr")->Fill(endVertex.Perp());
    H("EVz")->Fill(endVertex.Z());
    H("PEta")->Fill(p.Eta());
  }

  // Checking the lifetime
  Double_t fd = diffVertex.Mag();
  const Double_t c = TMath::C() * 1e-6;// * 1e3 / 1e9;//We need mm/ns in LHCb Units
  Double_t ltime = fd * (*D_MM) / (p.Mag() * c);
  H("Lifetime")->Fill(ltime / *D_BPVLTIME);
  
  // Keeping a histogram of acceptance ratio
  H("AcceptanceRatio")->Fill(acceptanceRatio);

  // Lifetime acceptance plot
  FillAcceptance(H("Acceptance"), acceptance);
  FillAcceptance(H("AcceptanceZ"), acceptance);

  // Checking the difference if we use the vertices direction instead
  auto  zAtCutFollowingVertices = GetZforRadius(radiusCut, primaryVertex, diffVertex);
  Double_t acceptanceV = *D_BPVLTIME * (zAtCutFollowingVertices - primaryVertex.Z()) / (endVertex.Z() - primaryVertex.Z());
  FillAcceptance(H("AcceptanceV"), acceptanceV);
  FillAcceptance(H("AcceptanceVZ"), acceptanceV);

  // And we're done...
  return kTRUE;
}

/**
 * @See header
 */
void RadialSelector::FillAcceptance(TH1D *hist, Double_t value) 
{
   
  TAxis *xaxis = hist->GetXaxis();
  Int_t cutoffbin = xaxis->FindBin(value);
  for (Int_t i = xaxis->GetFirst(); i <= cutoffbin; i++)
  {
      hist->AddBinContent(i);
  }
  // Set the number of entries as it isn't set by AddBinContent
  // If unset that creates problems for the merging in PROOF
  hist->SetEntries(hist->GetEntries() + 1);
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
  std::cout << "Hist event count  : " 
            << H("Count")->GetBinContent(H("Count")->GetXaxis()->FindBin("TOTAL")) 
            << std::endl;
  std::cout << "Tree total        : " 
            << H("Count")->GetBinContent(H("Count")->GetXaxis()->FindBin("TREE_TOTAL")) 
            << std::endl;
  std::cout << "Negative lifetime : " 
            << H("Count")->GetBinContent(H("Count")->GetXaxis()->FindBin("NEG_TIME")) 
            << std::endl;
  std::cout << "Negative ratio     : " 
            << H("Count")->GetBinContent(H("Count")->GetXaxis()->FindBin("NEG_ACCRATIO")) 
            << std::endl;
  std::cout << "Both negative      : " 
            << H("Count")->GetBinContent(H("Count")->GetXaxis()->FindBin("NEG_BOTH")) 
            << std::endl;
  std::cout << "EVz < PVz          : " 
            << H("Count")->GetBinContent(H("Count")->GetXaxis()->FindBin("EV_BEF_PV")) 
            << std::endl;
  std::cout << "EVz < PVz neg lifetime: " 
            << H("Count")->GetBinContent(H("Count")->GetXaxis()->FindBin("NEG_TIME_EV_BEF_PV")) 
            << std::endl;

  
  Int_t nbevents = H("Count")->GetBinContent(H("Count")->GetXaxis()->FindBin("TOTAL"));
  
  TCanvas *c1 = new TCanvas("c1","Histogram checks",200,10,700,900);
  c1->Divide(2,3);
  c1->cd(1);
  auto histMM = H("MM");
  histMM->Fit("gaus");
  histMM->Draw();
  c1->cd(3);
  auto histPVz = H("PVz");
  histPVz->Fit("gaus");
  histPVz->Draw();
  c1->cd(4);
  H("PVr")->Draw();
  c1->cd(5);
  auto histEVz = H("EVz");
  histEVz->SetStats(0);
  histEVz->DrawNormalized();
  auto histEVz_FromB = H("EVz_FromB");
  histEVz_FromB->SetStats(0);
  histEVz_FromB->SetLineColor(kGreen);
  histEVz_FromB->DrawNormalized("SAME");
  auto legendEVz = new TLegend(0.5,0.8,0.89,0.89);
  legendEVz->SetFillColor(0);
  legendEVz->AddEntry(histEVz);
  legendEVz->AddEntry(histEVz_FromB);
  legendEVz->Draw("SAME");

  c1->cd(6);
  auto histEVr = H("EVr");
  histEVr->SetStats(0);
  histEVr->DrawNormalized();
  auto histEVr_FromB = H("EVr_FromB");
  histEVr_FromB->SetStats(0);
  histEVr_FromB->SetLineColor(kGreen);
  histEVr_FromB->DrawNormalized("SAME");
  auto legendEVr = new TLegend(0.5,0.8,0.89,0.89);
  legendEVr->SetFillColor(0);
  legendEVr->AddEntry(histEVr);
  legendEVr->AddEntry(histEVr_FromB);
  legendEVr->Draw("SAME");

  c1->cd(2);
  auto histPEta = H("PEta");
  histPEta->SetStats(0);
  histPEta->DrawNormalized();
  auto histPEta_FromB = H("PEta_FromB");
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
  H("Lifetime")->DrawNormalized();
  c2->Update();
  
  TCanvas *c4 = new TCanvas("c4","Ratio between lifetime and max lifetime",200,10,700,900);
  H("AcceptanceRatio")->DrawNormalized();
  c4->Update();

  TCanvas *c3 = new TCanvas("c3","Acceptance",200,10,700,900);
  c3->SetLogy();
  c3->SetTitle("Lifetime acceptance (ns)");
  auto histAcceptance = H("Acceptance");
  histAcceptance->Scale(1.0/nbevents);
  histAcceptance->SetTitle("Lifetime acceptance effect due to 4mm radial cut");
  histAcceptance->SetStats(0);
  histAcceptance->GetXaxis()->SetTitle("time (ns)");
  histAcceptance->Draw();
  auto histAcceptanceV = H("AcceptanceV");
  histAcceptanceV->Scale(1.0/nbevents);
  histAcceptanceV->SetStats(0);
  histAcceptanceV->SetLineColor(kGreen);
  histAcceptanceV->Draw("SAME");
  auto legendLf = new TLegend(0.3,0.85,0.89,0.89);
  legendLf->SetFillColor(0);
  legendLf->AddEntry(histAcceptance);
  legendLf->AddEntry(histAcceptanceV);
  legendLf->Draw("SAME");
  c3->Update();

  TCanvas *c5 = new TCanvas("c5","Acceptance",200,10,700,900);
  c5->SetTitle("Lifetime acceptance (ns)");
  auto histAcceptanceZ = H("AcceptanceZ");
  histAcceptanceZ->Scale(1.0/nbevents);
  histAcceptanceZ->SetTitle("Lifetime acceptance effect due to 4mm radial cut");
  histAcceptanceZ->SetStats(0);
  histAcceptanceZ->GetXaxis()->SetTitle("time (ns)");
  histAcceptanceZ->Draw();
  auto histAcceptanceVZ = H("AcceptanceVZ");
  histAcceptanceVZ->Scale(1.0/nbevents);
  histAcceptanceVZ->SetStats(0);
  histAcceptanceVZ->SetLineColor(kGreen);
  histAcceptanceVZ->Draw("SAME");
  auto legendLfZ = new TLegend(0.3,0.85,0.89,0.89);
  legendLfZ->SetFillColor(0);
  legendLfZ->AddEntry(histAcceptanceZ);
  legendLfZ->AddEntry(histAcceptanceVZ);
  legendLfZ->Draw("SAME");
  c5->Update();
}
