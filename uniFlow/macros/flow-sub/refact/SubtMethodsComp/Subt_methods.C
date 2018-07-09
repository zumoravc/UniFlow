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

TH1D* PrepScalingFactor(TProfile* pMultRaw,TProfile* pMultBase);
Double_t GetScalingFactor(TProfile* pMultRaw, Int_t iCentRaw,TProfile* pMultBase, Int_t iCentBase);
TH1D* SubtractDn(TH1D* raw, TH1D* base, Double_t factor);
TH1D* SubtractCn(TH1D* raw, TH1D* base, TH1D* hFactor);
TH1D* MakeVn(TH1D* hDn, TH1D* hCn, Int_t iCent);

TH1D* Scale(TH1D* base, Double_t factor);
TH1D* ScaleCn(TH1D* base, TH1D* mult);


Color_t colBase = kBlue;
Color_t colRaw = kRed;
Color_t colSubt = kGreen+2;
// Color_t colSubt = colRaw;

Int_t markBase = kOpenCircle;
Int_t markBaseScaled = kFullCircle;
Int_t markRaw = kOpenSquare;
Int_t markSubt = kFullSquare;

const Int_t iNumCent = 4;
TString sCentLabel[iNumCent] = {"0-20%", "20-40%", "40-60%", "60-100%"};

TString sSpecies = "Charged";
TString sMethod = "GF_eventweighted";
TString sGap = "gap04";
Int_t iCentBase = 0; // centrality index; if -1 use same centrality as in Raw case

TString sPathRaw = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/pPb-run3-gaps-04-06-10-12/output_vn/%s/%s/",sGap.Data(),sMethod.Data());
TString sPathBase = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/pp-run3-gaps-04-06-10-12/output_vn_int/%s/%s/",sGap.Data(),sMethod.Data());
TString sPathOut = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/output-test/%s/%s/",sGap.Data(),sMethod.Data());

TString sFileIn = "Processed.root";
TString sFileOut = "Subtracted.root";
TString sFileMult = "Mult.root";
// ==================================================================================================================

void Subt_methods()
{
  // NOTE: conversion from cent index (0,..,iNumCent-1) to bin index (1,..,iNumCent)
  iCentBase++;

  // === LOADING INPUT ===
  // output files
  gSystem->mkdir(sPathOut,kTRUE);  // Making output folder
  TFile* fileOut = OpenFile(sPathOut+sFileOut,"RECREATE"); if(!fileOut) { return; }
  // multiplicities
  TFile* fileMultRaw = OpenFile(sPathRaw+"../"+sFileMult); if(!fileMultRaw) { return; }
  TFile* fileMultBase = OpenFile(sPathBase+"../"+sFileMult); if(!fileMultBase) { return; }
  // input files
  TFile* fileInRaw = OpenFile(sPathRaw+sFileIn); if(!fileInRaw) { return; }
  TFile* fileInBase = OpenFile(sPathBase+sFileIn); if(!fileInBase) { return; }

  printf("Files opened\n");

  // mult
  TProfile* pMult_Raw = (TProfile*) fileMultRaw->Get("fpRefsMult_rebin"); if(!pMult_Raw) { printf("ERROR pMultRaw not found!\n"); fileMultRaw->ls(); return; }
  pMult_Raw->SetName("pMult_Raw");
  TProfile* pMult_Base = (TProfile*) fileMultBase->Get("fpRefsMult_rebin"); if(!pMult_Base) { printf("ERROR pMultBase not found!\n"); fileMultBase->ls(); return; }
  pMult_Base->SetName("pMult_Base");
  TH1D* hMult_Scaled = PrepScalingFactor(pMult_Raw, pMult_Base);

  // // cn
  TH1D* hCn_Raw = LoadHisto(Form("hCum2_Refs_harm2_%s",sGap.Data()),fileInRaw); if(!hCn_Raw) { return; }
  TH1D* hCn_Base = LoadHisto(Form("hCum2_Refs_harm2_%s",sGap.Data()),fileInBase); if(!hCn_Base) { return; }
  TH1D* hCn_Sub = SubtractCn(hCn_Raw, hCn_Base, hFactor); if(!hCn_Sub) { return; }


  fileOut->cd();
  pMult_Raw->Write("pMult_Raw");
  pMult_Base->Write("pMult_Base");
  hFactor->Write();
  hCn_Raw->Write();
  hCn_Sub->Write();


  // dn  // one cent only now
  for(Int_t iCent(1); iCent < hCn_Raw->GetNbinsX()+1; ++iCent)
  {
    Int_t iBase = iCent;
    if(iCentBase > 0) { iBase = iCentBase; }

    TH1D* hDn_Raw = LoadHisto(Form("hCum2_%s_harm2_%s_cent%d",sSpecies.Data(),sGap.Data(),iCent-1),fileInRaw); if(!hCn_Raw) { return; }
    TH1D* hDn_Base = LoadHisto(Form("hCum2_%s_harm2_%s_cent%d",sSpecies.Data(),sGap.Data(),iCentBase-1),fileInBase); if(!hCn_Base) { return; }

    // === Getting scaling factor ===
    Double_t dFactor = hFactor->GetBinContent(iCent);

    // === SUBTRACTING ===
    TH1D* hDn_Sub = SubtractDn(hDn_Raw, hDn_Base, dFactor); if(!hDn_Sub) { return; }

    // === MAKING vn's ===
    TH1D* hVn_Raw = MakeVn(hDn_Raw, hCn_Raw, iCent); if(!hVn_Raw) { return; }
    printf("Raw vn done\n");
    TH1D* hVn_Sub = MakeVn(hDn_Sub, hCn_Sub, iCent); if(!hVn_Sub) { return; }
    printf("Sub vn done\n");

    fileOut->cd();
    hDn_Raw->Write();
    hVn_Raw->Write();
    hDn_Sub->Write();
    hVn_Sub->Write();
  }

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
TH1D* PrepScalingFactor(TProfile* pMultRaw, TProfile* pMultBase)
{
  if(!pMultRaw) { printf("ERROR-PrepScalingFactor: 'pMultRaw' does not found.\n"); return 0x0; }
  if(!pMultBase) { printf("ERROR-PrepScalingFactor: 'pMultBase' does not found.\n"); return 0x0; }

  TH1D* hFactor = (TH1D*) pMultRaw->ProjectionX("hMult_Scaled");
  hFactor->Reset();

  Double_t dContBase = pMultBase->GetBinContent(iCentBase);
  Double_t dErrBase = pMultBase->GetBinContent(iCentBase);

  for(Int_t bin(1); bin < hFactor->GetNbinsX()+1; ++bin)
  {
    if(iCentBase < 1) { dContBase = pMultBase->GetBinContent(bin); dErrBase = pMultBase->GetBinError(bin); }

    Double_t dContRaw = pMultRaw->GetBinContent(bin);
    Double_t dErrRaw = pMultRaw->GetBinError(bin);

    if(dContRaw <= 0.0) { hFactor->SetBinContent(bin,0.0); hFactor->SetBinError(bin,999.9);  continue; }

    Double_t dFactor = dContBase / dContRaw;
    Double_t dErrSq = TMath::Power(dErrBase/dContRaw,2.0) + TMath::Power( dContBase * dErrRaw / (dContRaw*dContRaw),2.0);

    hFactor->SetBinContent(bin,dFactor);
    hFactor->SetBinError(bin,TMath::Sqrt(dErrSq));
  }

  return hFactor;
}
// ==================================================================================================================
Double_t GetScalingFactor(TProfile* pMultRaw, Int_t iCentRaw,TProfile* pMultBase, Int_t iCentBase)
{
  if(!pMultRaw) { printf("ERROR-GetScalingFactor: 'pMultRaw' does not found.\n"); return -999.9; }
  if(!pMultBase) { printf("ERROR-GetScalingFactor: 'pMultBase' does not found.\n"); return -999.9; }

  if(iCentRaw < 0 || iCentRaw > pMultRaw->GetNbinsX()) { printf("ERROR-GetScalingFactor: Invalid 'iCentRaw'.\n"); return -999.9; }
  if(iCentBase < 0 || iCentBase > pMultBase->GetNbinsX()) { printf("ERROR-GetScalingFactor: Invalid 'iCentBase'.\n"); return -999.9; }

  Double_t dRaw = pMultRaw->GetBinContent(iCentRaw);
  Double_t dBase = pMultBase->GetBinContent(iCentBase);

  Double_t dFactor = dBase / dRaw;

  return dFactor;
}
// ==================================================================================================================
TH1D* SubtractCn(TH1D* raw, TH1D* base, TH1D* hFactor)
{
  if(!raw) { printf("ERROR-SubtractCn: Hist 'raw' does not found.\n"); return 0x0; }
  if(!base) { printf("ERROR-SubtractCn: Hist 'base' does not found.\n"); return 0x0; }
  if(!hFactor) { printf("ERROR-SubtractCn: Hist 'hFactor' does not found.\n"); return 0x0; }

  TH1D* sub = (TH1D*) raw->Clone(Form("%s_sub",raw->GetName()));
  if(!sub) { printf("ERROR-SubtractCn: Hist 'sub' does not cloned properly.\n"); return 0x0; }

  for(Int_t bin(1); bin < sub->GetNbinsX()+1; ++bin)
  {
    Double_t con_raw = raw->GetBinContent(bin);
    Double_t err_raw = raw->GetBinError(bin);
    Double_t con_base = base->GetBinContent(iCentBase);
    Double_t err_base = base->GetBinError(iCentBase);

    Double_t factor = hFactor->GetBinContent(bin);
    Double_t factorErr = hFactor->GetBinError(bin);

    Double_t con = con_raw - factor * con_base;
    Double_t errSq = TMath::Power(err_raw,2.0) + TMath::Power(err_base*factor,2.0) + TMath::Power(con_base*factorErr,2.0);

    sub->SetBinContent(bin, con);
    sub->GetBinError(bin, TMath::Sqrt(errSq));
  }

  return sub;
}
// ==================================================================================================================
TH1D* SubtractDn(TH1D* raw, TH1D* base, Double_t factor)
{
  if(!raw) { printf("ERROR-SubtractDn: Hist 'raw' does not found.\n"); return 0x0; }
  if(!base) { printf("ERROR-SubtractDn: Hist 'base' does not found.\n"); return 0x0; }

  TH1D* sub = (TH1D*) base->Clone(Form("%s_sub",raw->GetName()));
  if(!sub) { printf("ERROR-SubtractDn: Hist 'sub' does not cloned properly.\n"); return 0x0; }

  for(Int_t bin(1); bin < sub->GetNbinsX()+1; ++bin)
  {
    Double_t con_raw = raw->GetBinContent(bin);
    Double_t err_raw = raw->GetBinError(bin);
    Double_t con_base = base->GetBinContent(bin);
    Double_t err_base = base->GetBinError(bin);

    sub->SetBinContent(bin, con_raw - factor * con_base);
    sub->GetBinError(bin, TMath::Sqrt(err_raw*err_raw + factor*factor*err_base*err_base));
  }

  return sub;
}
// ==================================================================================================================
TH1D* MakeVn(TH1D* hDn, TH1D* hCn, Int_t iCent)
{
  if(!hDn) { printf("ERROR-MakeVn: Hist 'hDn' does not found.\n"); return 0x0; }
  if(!hCn) { printf("ERROR-MakeVn: Hist 'hCn' does not found.\n"); return 0x0; }

  TH1D* hVn = (TH1D*) hDn->Clone(Form("%s_vn",hDn->GetName()));
  if(!hVn) { printf("ERROR-MakeVn: Hist 'hVn' does not cloned properly.\n"); return 0x0; }
  hVn->Reset();

  Double_t dConCn = hCn->GetBinContent(iCent);
  Double_t dErrCn = hCn->GetBinError(iCent);
  if(dConCn <= 0.0) { printf("ERROR-MakeVn: Division by zero (cent %d : %f)\n",iCent, dConCn); hCn->Draw(); return hVn; }


  for(Int_t bin(0); bin < hDn->GetNbinsX()+1; ++bin)
  {
    Double_t dConDn = hDn->GetBinContent(bin);
    Double_t dErrDn = hDn->GetBinError(bin);

    Double_t dConVn = dConDn / TMath::Sqrt(dConCn);
    Double_t dErrVnSq = dErrDn*dErrDn / dConCn + 0.25*dConDn*dConDn*dErrCn*dErrCn/(dConCn*dConCn*dConCn);

    hVn->SetBinContent(bin, dConVn);
    hVn->SetBinError(bin, TMath::Sqrt(dErrVnSq));
  }

  return hVn;
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
TH1D* Scale(TH1D* base, Double_t factor)
{
  if(!base) { printf("ERROR-Scale: Hist 'base' does not found.\n"); return 0x0; }

  TH1D* scale = (TH1D*) base->Clone(Form("%s_scaled",base->GetName()));
  if(!scale) { printf("ERROR-scale: Hist 'scale' does not cloned properly.\n"); return 0x0; }

  for(Int_t bin(1); bin < base->GetNbinsX()+1; ++bin)
  {
    scale->SetBinContent(bin, factor * base->GetBinContent(bin));
    scale->SetBinError(bin, TMath::Abs(factor * base->GetBinError(bin)));
  }

  return scale;
}
// ==================================================================================================================
TH1D* ScaleCn(TH1D* base, TH1D* mult)
{
  if(!base) { printf("ERROR-ScaleCn: Hist 'base' does not found.\n"); return 0x0; }
  if(!mult) { printf("ERROR-ScaleCn: Hist 'mult' does not found.\n"); return 0x0; }

  TH1D* scaled = (TH1D*) base->Clone(Form("%s_scaled",base->GetName()));

  for(Int_t bin(1); bin < base->GetNbinsX()+1; ++bin)
  {
    Double_t con =  base->GetBinContent(bin);
    Double_t err =  base->GetBinError(bin);
    Double_t mult_con =  mult->GetBinContent(bin);
    Double_t mult_err =  mult->GetBinError(bin);

    scaled->SetBinContent(bin, con * mult_con*mult_con );
    scaled->SetBinError(bin, TMath::Power( mult_con*mult_con*err, 2.0) + TMath::Power(2.0*con*mult_con*mult_err, 2.0));
  }

  return scaled;
}
// ==================================================================================================================
