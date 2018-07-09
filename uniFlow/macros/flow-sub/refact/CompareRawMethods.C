/* CompareRawMethods
 *
 * Compare raw (unsubtracted) result obtained from various flow-calculations methods.
 *
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2018
 */

#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

TFile* OpenFile(TString sFileName, TString sMode = "READ");
TH1D* LoadHisto(TString sHistName, TFile* file);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);

Color_t colors[] = {kGreen+2, kBlue, kRed, kBlack, kMagenta+1};

void CompareRawMethods()
{
  const Int_t iNumMethods = 3;
  TString sMethods[iNumMethods] = {"GF_eventweighted","GF_noneventweighted","SP_nonscaled_noneventweighted"};

  TString sSpecies = "Charged";
  TString sGap = "gap08";

  Double_t dMultBinning[] = {0,150};
  // Double_t dMultBinning[] = {20,25,30,35,40,45,50,55,60,65,70,75,80};
  const Int_t iNumMult = sizeof(dMultBinning)/sizeof(dMultBinning[0]);

  TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-NchRFP-gap08/output_vn_unit/";
  TString sOutputPath = sInputPath+"/CompareRawMethods/";
  // ##########################################################################################################
  gSystem->mkdir(sOutputPath.Data(),kTRUE);


  TLegend* leg = new TLegend(0.12,0.7,0.6,0.89);
  leg->SetHeader("p-Pb");
  leg->SetBorderSize(0.);
  leg->SetFillColor(0);

  // Plotting RFPs flow

  TCanvas* canRefs = new TCanvas("canRefs","canRefs",800,800);
  canRefs->cd();
  TH1* frame_canRefs = (TH1*) gPad->DrawFrame(dMultBinning[0],0.,dMultBinning[iNumMult-1],2.0);
  frame_canRefs->SetTitle(Form("Refs v_{2}{2,%s} (0-20%% V0A); N_{RFP}; v_{2}{2}",sGap.Data()));

  for(Int_t iMethod(0); iMethod < iNumMethods; ++iMethod)
  {
    TString sHistoName = Form("hFlow2_Refs_harm2_%s",sGap.Data());

    TFile* fileInput = OpenFile(sInputPath+sMethods[iMethod]+"/Processed.root","READ"); if(!fileInput) return;
    TH1D* flow = LoadHisto(sHistoName,fileInput); if(!flow) return;

    canRefs->cd();
    StyleHist(flow,colors[iMethod],kFullCircle);
    flow->Draw("same");
    leg->AddEntry(flow,sMethods[iMethod].Data(),"pel");
  }
  canRefs->cd();
  leg->Draw();
  canRefs->SaveAs(Form("%s/refs.pdf",sOutputPath.Data()),"pdf");

  // Plotting pt-diff flow

  for(Int_t iMult(0); iMult < iNumMult-1; ++iMult)
  {
    TString sMultLabel = Form("%2.0f-%2.0f",dMultBinning[iMult],dMultBinning[iMult+1]);

    TCanvas* canDiff = new TCanvas(Form("canDiff%d",iMult),Form("canDiff%d",iMult),800,800);
    canDiff->cd();
    TH1* frame_canDiff = (TH1*) gPad->DrawFrame(0.0,0.0,10.0,0.35);
    frame_canDiff->SetTitle(Form("%s v_{2}{2,%s} N_{RFP} %s (0-20%% V0A) ; p_{T} (GeV/c); v_{2}{2}",sSpecies.Data(),sGap.Data(),sMultLabel.Data()));

    for(Int_t iMethod(0); iMethod < iNumMethods; ++iMethod)
    {
      TString sHistoName = Form("hFlow2_%s_harm2_%s_cent%d",sSpecies.Data(),sGap.Data(),iMult);

      TFile* fileInput = OpenFile(sInputPath+sMethods[iMethod]+"/Processed.root","READ"); if(!fileInput) return;
      TH1D* flow = LoadHisto(sHistoName,fileInput); if(!flow) return;

      canDiff->cd();
      StyleHist(flow,colors[iMethod],kFullCircle);
      flow->Draw("same");
      // leg->AddEntry(flow,sMethods[iMethod].Data(),"pel");
    }
    canDiff->cd();
    leg->Draw();
    canDiff->SaveAs(Form("%s/%s_mult%d.pdf",sOutputPath.Data(),sSpecies.Data(),iMult),"pdf");

  }

  // printf("%s\n",sMultLabel.Data());

  return;
}
// ==================================================================================================================
TFile* OpenFile(TString sFileName, TString sMode)
{
  TFile* file = TFile::Open(sFileName.Data(),sMode.Data());
  if(!file) { printf("ERROR: Input file '%s' not found.\n",sFileName.Data()); return 0x0; }

  return file;
}
// ==================================================================================================================
TH1D* LoadHisto(TString sHistName, TFile* file)
{
  if(!file) { printf("ERROR-LoadHisto: File does not found.\n"); return 0x0; }

  TH1D* hist = (TH1D*) file->Get(sHistName.Data());
  if(!hist) { printf("ERROR-LoadHisto: Histo '%s' not found\n",sHistName.Data()); file->ls(); return 0x0; }

  return hist;
}
// ==================================================================================================================
void StyleHist(TH1* hist, Color_t color, Style_t markerStyle)
{
  if(!hist) { printf("ERROR-DrawHist: Hist does not found.\n"); return; }
  hist->SetLineColor(color);
  // hist->SetLineStyle(color);
  // hist->SetLineStyle(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  return;
};
// ==================================================================================================================
