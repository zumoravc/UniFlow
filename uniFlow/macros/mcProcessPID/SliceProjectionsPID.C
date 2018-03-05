// =================================================================================================
// PlotPIDProjections.C
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
//
// Load PID QA 3D histograms with nSigma TPC vs nSigma TOF vs pT
// and slice it up wrt. pt
//
// =================================================================================================

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMath.h"
#include "TColor.h"
#include "TLine.h"

void SliceProjectionsPID()
{
  TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/mc-pid/pPb-LHC17f2a_cent_woSDD_fix-run3";
  TString sOutputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/mc-pid/pPb-LHC17f2a_cent_woSDD_fix-run3";
  TString sTaskName = "UniFlow_3sigma";

  const Int_t iNumSpecies = 3;
  TString sSpeciesName[iNumSpecies] = {"Pion","Kaon","Proton"};
  Color_t colors[] = {kRed, kBlue+2, kGreen+1, kViolet-1};

  Double_t dPtBinsEdges[] = {0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,5.0,6.0,7.0,8.0,9.0,10.0,20.0};
  Int_t iNumPtBins = sizeof(dPtBinsEdges) / sizeof(dPtBinsEdges[0]);


  // ==================================================================================================================

  gSystem->mkdir(sOutputPath.Data(),kTRUE);
  TFile* fileOut = TFile::Open(Form("%s/outSlideProjectionsPID.root",sOutputPath.Data()), "RECREATE");
  if(!fileOut) { printf("ERROR: Input file not found!\n"); return; }

  TFile* fileIn = TFile::Open(Form("%s/AnalysisResults.root",sInputPath.Data()), "READ");
  if(!fileIn) { printf("ERROR: Input file not found!\n"); return; }

  fileIn->cd(sTaskName.Data());
  TList* listInPID = (TList*) gDirectory->Get(Form("QA_PID_%s",sTaskName.Data()));
  if(!listInPID) { printf("ERROR: Input 'listInPID' not found"); gDirectory->ls(); return; }

  // listInPID->ls();

  TCanvas* can = new TCanvas("can","can");

  for(Int_t iSpecies(0); iSpecies < iNumSpecies; ++iSpecies)
  {
    TH3D* hist3d = (TH3D*) listInPID->FindObject(Form("fh3QAPIDnSigmaTPCTOFPt%s_Before",sSpeciesName[iSpecies].Data()));
    if(!hist3d) { printf("ERROR: Input TH3D 'fh3QAPIDnSigmaTPCTOFPt%s_Before' not found", sSpeciesName[iSpecies].Data()); listInPID->ls(); return; }
    hist3d->Draw();

    for(Int_t iPt(0); iPt < iNumPtBins-1; ++iPt)
    {
      // lower bin
      Int_t iBinLow = hist3d->GetZaxis()->FindFixBin(dPtBinsEdges[iPt]);
       // upper bin (included)
      Int_t iBinHigh = hist3d->GetZaxis()->FindFixBin(dPtBinsEdges[iPt+1])-1;

      TH3D* hist3dTemp = (TH3D*) hist3d->Clone();
      hist3dTemp->GetZaxis()->SetRange(iBinLow,iBinHigh);
      Double_t dPtLow = hist3dTemp->GetZaxis()->GetBinLowEdge(iBinLow);
      Double_t dPtHigh = hist3dTemp->GetZaxis()->GetBinUpEdge(iBinHigh);

      printf("bin %d-%d : pt %g-%g\n", iBinLow,iBinHigh,dPtLow,dPtHigh);

      TH2D* h2_project = (TH2D*) hist3dTemp->Project3D("xy");
      h2_project->SetTitle(Form("nSigma TPC vs TOF (%s) %.1f < p_{T} < %.1f",sSpeciesName[iSpecies].Data(),dPtLow,dPtHigh));
      // h2_project->Draw("colz");

      fileOut->cd();
      h2_project->Write(Form("proj_%s_pt_%d_%d",sSpeciesName[iSpecies].Data(),iBinLow,iBinHigh));
    }
  }

  return;
}
