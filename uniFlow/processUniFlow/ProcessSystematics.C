/* ProcessSystematics
 *
 * Macro for processing results of systematical test.
 * Producing the ratios and estimating the overall uncertainty.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

const Short_t iNumFiles = 2;
const char* sTaskTag[iNumFiles] = {"vtx8","vtx9"};

const Short_t iNumSpecies = 1;
const char* sSpecies[iNumSpecies] = {"hFlow2_Charged_harm2_gap08_cent0"};

const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst/Vtx_z/syst";
const char* sOutputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst/Vtx_z/syst";

const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7};
const Int_t colors[]     = {kGreen+3, kBlue+1, kMagenta+1,kRed+1,kOrange-1,kCyan+2,kYellow+2};
const Int_t markers[]    = {kFullCircle, kFullSquare,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullCross,kFullDiamond,kFullStar,kOpenStar};

TCanvas* canSuper;
TCanvas* canRatio;
TCanvas* canBarlow;

void ProcessList(TList* list = 0x0);
TH1D* MakeRatio(TH1D* hBase = 0x0, TH1D* hSyst = 0x0, Bool_t bCorelated = kFALSE);
TH1D* DoBarlowTest(TH1D* hBase = 0x0, TH1D* hSyst = 0x0, Bool_t bCorelated = kFALSE);

void ProcessSystematics()
{
  canSuper = new TCanvas("Super","canSuper",800,800);
  canRatio = new TCanvas("canRatio","canRatio",800,800);
  canBarlow = new TCanvas("canBarlow","canBarlow",800,800);

  TList* list = new TList();

  TFile* fInputFile = 0x0;
  TH1D* hTemp = 0x0;

  // loop over files
  for(Short_t iFile(0); iFile < iNumFiles; iFile++)
  {
    fInputFile = new TFile(Form("%s/%s/UniFlow_%s.root",sInputPath,sTaskTag[iFile],sTaskTag[iFile]),"READ");
    if(!fInputFile) { printf("InputFile (%d) not found!\n",iFile); return; }

    fInputFile->ls();

    // loop over species (histos of interests)
    for(Short_t iSpecies(0); iSpecies < iNumSpecies; iSpecies++)
    {
      hTemp = (TH1D*) fInputFile->Get(sSpecies[iSpecies]);
      if(!hTemp) { printf("Histos not found: file %d | species %d\n",iFile,iSpecies); return; }
      list->Add(hTemp);
    }
  }

  printf("List entries: %d\n",list->GetEntries());
  ProcessList(list);

  return;
}
//_____________________________________________________________________________
void ProcessList(TList* list)
{
  if(!list) { printf("ProcessList::List not found\n"); return; }
  const Short_t iEntries = list->GetEntries();
  if(iEntries != iNumFiles) { printf("Different number of entries (%d) than expected (%d)\n",iEntries,iNumFiles); return; }

  canSuper->cd();

  TH1D* hBase = (TH1D*) list->At(0);
  if(!hBase) { printf("ProcessList::Baseline histo not found\n"); return; }

  hBase->SetLineColor(colors[0]);
  hBase->SetMarkerColor(colors[0]);
  hBase->SetMarkerStyle(markers[0]);
  hBase->SetFillColor(colors[0]);

  hBase->SetMinimum(0.);
  hBase->SetMaximum(1.);
  hBase->Draw();

  TH1D* hTemp = 0x0;
  TH1D* hRatio = 0x0;
  // TH1D* hRatio2 = 0x0;
  TH1D* hBarlow = 0x0;

  for(Short_t iFile(1); iFile < iNumFiles; iFile++)
  {
    hTemp = (TH1D*) list->At(iFile);
    if(!hTemp) { printf("ProcessList::Histo (%d) not found\n",iFile); return; }

    hTemp->SetLineColor(colors[iFile]);
    hTemp->SetMarkerColor(colors[iFile]);
    hTemp->SetMarkerStyle(markers[iFile]);
    hTemp->SetFillColor(colors[iFile]);

    canSuper->cd();
    hTemp->Draw("same");

    hRatio = MakeRatio(hBase,hTemp,1);
    // hRatio->SetLineColor(colors[iFile]);
    // hRatio->SetMarkerColor(colors[iFile]);
    // hRatio->SetMarkerStyle(markers[iFile]);
    // hRatio->SetFillColor(colors[iFile]);

    hBarlow = DoBarlowTest(hBase,hTemp,1);

    canRatio->cd();
    hRatio->Draw("hist p e1");

    canBarlow->cd();
    hBarlow->Draw("hist p");
  }

  // fitting the ratios
  canRatio->cd();
  TF1* fitConst = new TF1("fitConst","[0]",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
  fitConst->SetLineColor(kRed);
  hRatio->Fit("fitConst","R");
  printf("Fit: %g ± %g (Chi2 %g)\n",fitConst->GetParameter(0),fitConst->GetParError(0),fitConst->GetChisquare());

  TLine* lineUnity = new TLine();
  lineUnity->SetLineColor(kBlack);
  lineUnity->SetLineStyle(7);

  fitConst->Draw("same");
  lineUnity->DrawLine(hRatio->GetXaxis()->GetXmin(),1.,hRatio->GetXaxis()->GetXmax(),1.);

  canBarlow->cd();
  lineUnity->DrawLine(hBarlow->GetXaxis()->GetXmin(),1.,hBarlow->GetXaxis()->GetXmax(),1.);

  return;
}
//_____________________________________________________________________________
TH1D* MakeRatio(TH1D* hBase, TH1D* hSyst, Bool_t bCorelated)
{
  if(!hBase) { printf("MakeRatio::Baseline histo not found\n"); return 0x0; }
  if(!hSyst) { printf("MakeRatio::Systematics histo not found\n"); return 0x0; }

  TH1D* hRatio = (TH1D*) hSyst->Clone("hRatio");
  hRatio->Reset();


  // checking the binning
  const Short_t iNumBins = hRatio->GetNbinsX();
  const Short_t iNumBinsBase = hBase->GetNbinsX();
  const Short_t iNumBinsSyst = hSyst->GetNbinsX();

  if( iNumBins != iNumBinsBase || iNumBins != iNumBinsSyst ) { printf("MakeRatio::Different number of bins\n"); return 0x0; }

  // calculating ratio
  Double_t dContBase = 0, dErrBase = 0;
  Double_t dContSyst = 0, dErrSyst = 0;
  Double_t dContRatio = 0, dErrRatio = 0;
  for(Short_t iBin(1); iBin < iNumBins+1; iBin++)
  {
    dContBase = hBase->GetBinContent(iBin);
    dErrBase = hBase->GetBinError(iBin);
    dContSyst = hSyst->GetBinContent(iBin);
    dErrSyst = hSyst->GetBinError(iBin);

    dContRatio = dContBase / dContSyst;
    dErrRatio = TMath::Power(dErrBase/dContSyst, 2) + TMath::Power( dErrSyst*dContBase/(dContSyst*dContSyst), 2);
    // printf("Err (before) : %g | ", TMath::Sqrt(dErrRatio));

    if(bCorelated) dErrRatio += (2*dContBase*dErrBase*dErrSyst/TMath::Power(dContSyst,3));
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    hRatio->SetBinContent(iBin,dContRatio);
    hRatio->SetBinError(iBin,TMath::Sqrt(dErrRatio));
  }

  hRatio->SetMinimum(0.);
  hRatio->SetMaximum(2.);

  return hRatio;
}
//_____________________________________________________________________________
TH1D* DoBarlowTest(TH1D* hBase, TH1D* hSyst, Bool_t bCorelated)
{
  if(!hBase) { printf("DoBarlowTest::Baseline histo not found\n"); return 0x0; }
  if(!hSyst) { printf("DoBarlowTest::Systematics histo not found\n"); return 0x0; }

  TH1D* hBarlow = (TH1D*) hSyst->Clone("hBarlow");
  hBarlow->Reset();

  // checking the binning
  const Short_t iNumBins = hBarlow->GetNbinsX();
  const Short_t iNumBinsBase = hBase->GetNbinsX();
  const Short_t iNumBinsSyst = hSyst->GetNbinsX();

  if( iNumBins != iNumBinsBase || iNumBins != iNumBinsSyst ) { printf("DoBarlowTest::Different number of bins\n"); return 0x0; }

  // calculating ratio
  Double_t dContBase = 0, dErrBase = 0;
  Double_t dContSyst = 0, dErrSyst = 0;
  Double_t dBarlow = 0, dDenom = 0;
  for(Short_t iBin(1); iBin < iNumBins+1; iBin++)
  {
    dContBase = hBase->GetBinContent(iBin);
    dErrBase = hBase->GetBinError(iBin);
    dContSyst = hSyst->GetBinContent(iBin);
    dErrSyst = hSyst->GetBinError(iBin);

    if(bCorelated) { dDenom = TMath::Sqrt(TMath::Abs(TMath::Power(dErrBase,2) - TMath::Power(dErrSyst,2))); }
    else { dDenom = TMath::Sqrt(TMath::Abs(TMath::Power(dErrBase,2) + TMath::Power(dErrSyst,2))); }

    dBarlow = (dContBase - dContSyst) / dDenom;

    hBarlow->SetBinContent(iBin,dBarlow);
    hBarlow->SetBinError(iBin,0);
  }

  // hBarlow->SetMinimum(0.);
  // hBarlow->SetMaximum(10.);

  return hBarlow;
}
