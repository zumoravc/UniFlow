/*
* Compare effect of NUA (difference between applying and not appliying NUA) on cumulant level.
* Taken from Processed.root files.
*
* author : Vojtech Pacik (NBI) vojtech.pacik@cern.ch
*/

#include "TFile.h"
#include "TList.h"
#include "TH1.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TString.h"
#include "TCanvas.h"
#include "TMath.h"

void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle, Bool_t showStats = kFALSE);
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);

Bool_t bProcessSingleCase(TString sFileWithNUA, TString sFileNoNUA, TString sOutputPath, TString sGap, TString sHistName);

Color_t colWithNUA = kGreen+2;
Color_t colNoNUA = kRed;
Int_t markWithNUA = kFullCircle;
Int_t markNoNUA = kOpenSquare;

void CompareNUAeffect()
{
  TString sGap = "gap04";
  TString sSpecies[] = {"Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};
  TString sFileWithNUA = Form("/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/output_binning-3-test/%s/Processed.root",sGap.Data());
  TString sFileNoNUA = Form("/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt/output_binning-3-test/%s/Processed.root",sGap.Data());
  // TString sHistName = Form("hCum2_Refs_harm2_%s",sGap.Data());
  TString sOutputPath = Form("/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/outcome/NUAeffect/merged-pPb-16qt/");


  for(Int_t iSpecies(0); iSpecies < (sizeof(sSpecies)/sizeof(sSpecies[0])); ++iSpecies)
  {
    TString sHistName = Form("hCum2_%s_harm2_%s",sSpecies[iSpecies].Data(),sGap.Data());

    for(Int_t iCent(0); iCent < 4; ++iCent)
    {
      TString sHistNameCent = Form("%s_cent%d", sHistName.Data(), iCent);
      if(sSpecies[iSpecies].EqualTo("Refs")) { sHistNameCent = sHistName; }

      Bool_t bProc = bProcessSingleCase(sFileWithNUA, sFileNoNUA, sOutputPath, sGap, sHistNameCent);
      if(!bProc) { printf("ERROR: ProcessSingleCase (iCent %d) failed!\n",iCent); return; }
    }
  }

  return;
}
// ==================================================================================================================
Bool_t bProcessSingleCase(TString sFileWithNUA, TString sFileNoNUA, TString sOutputPath, TString sGap, TString sHistName)
{
  // ========================================
  gSystem->mkdir(sOutputPath.Data(),kTRUE);

  TFile* fileWithNUA = TFile::Open(sFileWithNUA.Data(),"READ"); if(!fileWithNUA) { printf("ERROR: no fileWithNUA\n"); return kFALSE; }
  TFile* fileNoNUA = TFile::Open(sFileNoNUA.Data(),"READ"); if(!fileNoNUA) { printf("ERROR: no fileNoNUA\n"); return kFALSE; }

  TH1D* histWithNUA = (TH1D*) fileWithNUA->Get(sHistName.Data()); if(!histWithNUA) { printf("ERROR: no histWithNUA '%s'\n",sHistName.Data()); fileWithNUA->ls(); return kFALSE; }
  TH1D* histNoNUA = (TH1D*) fileNoNUA->Get(sHistName.Data()); if(!histNoNUA) { printf("ERROR: no histNohNUA '%s'\n",sHistName.Data()); fileNoNUA->ls(); return kFALSE; }

  TH1D* ratio = DivideHistos(histNoNUA,histWithNUA,kTRUE);

  StyleHist(histWithNUA,colWithNUA,markWithNUA);
  StyleHist(histNoNUA,colNoNUA,markNoNUA);

  TCanvas* can = new TCanvas("can","can",600,600);
  can->Divide(2,1);
  can->cd(1);
  histWithNUA->DrawCopy();
  histNoNUA->DrawCopy("same");
  can->cd(2);
  ratio->DrawCopy();

  can->SaveAs(Form("%s/%s.pdf",sOutputPath.Data(),sHistName.Data()),"pdf");

  return kTRUE;
}
// ==================================================================================================================
void StyleHist(TH1* hist, Color_t color, Style_t markerStyle, Bool_t showStats)
{
  if(!hist) { printf("ERROR-DrawHist: Hist does not found.\n"); return; }
  hist->SetStats(showStats);
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerStyle);
  return;
};
// ==================================================================================================================
TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor)
{
  if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

  Int_t binsNom = nom->GetNbinsX();
  Int_t binsDenom = denom->GetNbinsX();

  // if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  TH1D* ratio = (TH1D*) nom->Clone(Form("Ratio_%s_%s",nom->GetName(),denom->GetName()));
  ratio->Reset();

  Double_t dContNom = 0, dErrNom = 0;
  Double_t dContDenom = 0, dErrDenom = 0;
  Double_t dContRatio = 0, dErrRatio = 0;
  for(Short_t iBin(1); iBin < binsDenom+1; iBin++)
  {
    if(iBin > binsNom) break;

    dContNom = nom->GetBinContent(iBin);
    dErrNom = nom->GetBinError(iBin);
    dContDenom = denom->GetBinContent(iBin);
    dErrDenom = denom->GetBinError(iBin);

    if(dContDenom == 0.0) continue;

    dContRatio =  dContNom / dContDenom;
    dErrRatio = TMath::Power(dErrNom/dContDenom, 2) + TMath::Power( dErrDenom*dContNom/(dContDenom*dContDenom), 2);
    // printf("Err (before) : %g | ", TMath::Sqrt(dErrRatio));

    if(bCor) dErrRatio -= (2*dContNom*dErrDenom*dErrNom/TMath::Power(dContDenom,3));
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    ratio->SetBinContent(iBin,dContRatio);
    ratio->SetBinError(iBin,TMath::Sqrt(dErrRatio));
  }

  return ratio;
}
// ==================================================================================================================
