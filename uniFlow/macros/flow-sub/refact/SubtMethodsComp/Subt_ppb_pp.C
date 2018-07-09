#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

void Subt_ppb_pp_species_gap(TString sGap="gap08");

TFile* OpenFile(TString sFileName, TString sMode = "READ");
TH1D* LoadHisto(TString sHistName, TFile* file);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);

TH1D* Scale(TH1D* base, Double_t factor, Double_t factorErrSq);
TH1D* ScaleCn(TH1D* base, TH1D* mult);
TH1D* Subtract(TH1D* raw, TH1D* base);
TH1D* Subtract_new_cn(TH1D* raw, TH1D* base, TProfile* mult_raw, TProfile* mult_base);
TH1D* Subtract_new_dn(TH1D* raw, TH1D* base,  Double_t dFactor, Double_t dFactor_err_sq);


// colors for centrality

Color_t colBase = kBlue;
Color_t colRaw = kRed;
Color_t colSubt = kGreen+2;
// Color_t colSubt = colRaw;

Int_t markBase = kOpenCircle;
Int_t markBaseScaled = kFullCircle;
Int_t markRaw = kOpenSquare;
Int_t markSubt = kFullSquare;

// TString sSpecies_list[] = {"Charged"};
// TString sSpecies_list[] = {"Charged","Pion","Kaon","Proton","K0s"};
// TString sSpecies_list[] = {"K0s","Kaon"};
TString sSpecies_list[] = {"Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};
Int_t iNumSpecies = sizeof(sSpecies_list) / sizeof(sSpecies_list[0]);

TFile* fileOutSubt = 0x0;
TString sOutPath;

void Subt_ppb_pp()
{
  TString sGap[] = {"gap00","gap04","gap08"};
  // TString sGap[] = {"gap08"};

  Int_t iNumGaps = sizeof(sGap) / sizeof(sGap[0]);

  for(Int_t g(0); g < iNumGaps; ++g)
  {
    // sOutPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/test/output_binning/"+sGap[g]+"/";
    sOutPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt/subt/output_binning-3/"+sGap[g]+"/";
    gSystem->mkdir(sOutPath.Data(),kTRUE);
    fileOutSubt = TFile::Open(Form("%sSubtracted.root",sOutPath.Data()),"RECREATE");

    Subt_ppb_pp_species_gap(sGap[g]);
    // Subt_ppb_pp_species_gap(sSpecies_list[i], sGap[g]);
  }

  return;
}


void Subt_ppb_pp_species_gap(TString sGap)
{
  // TString sMethod = "GF_eventweighted";
  // TString sOutputTag = "output_vn";
  // TString sOutputTagInt = sOutputTag + "_int";


  TString sGapBase = sGap;
  // TString sGapBase = "gap00";
  TString sGapRaw = sGapBase;

  TString sInFileRaw = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt/output_binning-3-test/" + sGapBase; // + "/" + sMethod;
  TString sInFileBaseInt = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pp-16kl-nua/output_binning-3/" + sGapBase; // + "/"; + sMethod;
  // TString sOutFolder = sInFileRaw+"/"+sMethod+"/pPb_pp_subt_"+sGapRaw+"/"+sSpecies;

  TString sOutFolder = sOutPath+"/test/";
  // TString sOutFolder = sOutPath+sSpecies;
  TString sOutFile = sOutFolder+"/Subt_results.root";

  const Int_t iNumCent = 5;
  TString sCentLabel[iNumCent] = {"0-10%","10-20%", "20-40%", "40-60%", "60-100%"};

  // ==================================================================================================================

  // === LOADING INPUT ===
  // multiplicities
  // TFile* fileInRaw_Mult = OpenFile(sInFileRaw+"/Mult.root"); if(!fileInRaw_Mult) { return; }
  // TProfile* hRaw_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInRaw_Mult); if(!hRaw_Mult) { return; }

  // TFile* fileInBase_MultInt = OpenFile(sInFileBaseInt+"/Mult.root"); if(!fileInBase_MultInt) { return; }
  // TProfile* hBase_MultInt = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBase_MultInt); if(!hBase_MultInt) { return; }

  // output files
  gSystem->mkdir(sOutFolder,kTRUE);  // Making output folder
  TFile* fileOut = OpenFile(sOutFile,"RECREATE"); if(!fileOut) { return; }

  // input files
  TFile* fileInRaw = OpenFile(sInFileRaw+"/Processed.root"); if(!fileInRaw) { return; }
  TFile* fileInBaseInt = OpenFile(sInFileBaseInt+"/Processed.root"); if(!fileInBaseInt) { return; }

  // mult
  TProfile* hRaw_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInRaw); if(!hRaw_Mult) { return; }
  TProfile* hBase_MultInt = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBaseInt); if(!hBase_MultInt) { return; }

  TH1D* hScaled_Mult = (TH1D*) hRaw_Mult->ProjectionX("hScaled_Mult");
  hScaled_Mult->Reset();

  StyleHist(hRaw_Mult,colRaw,markRaw);
  StyleHist(hBase_MultInt,colBase,markBase);
  StyleHist(hScaled_Mult,colSubt,markSubt);

  Double_t dMultBase = hBase_MultInt->GetBinContent(1);
  Double_t dMultBaseErr = hBase_MultInt->GetBinError(1);

  for(Int_t bin(1); bin < hScaled_Mult->GetNbinsX()+1; ++bin)
  {
    Double_t dMultRaw = hRaw_Mult->GetBinContent(bin);
    Double_t dMultRawErr = hRaw_Mult->GetBinError(bin);

    if(dMultRaw <= 0.0) { hScaled_Mult->SetBinContent(bin,0.0); hScaled_Mult->SetBinError(bin,999.9);  continue; }

    Double_t dErrSq = TMath::Power(dMultBaseErr/dMultRaw,2.0) + TMath::Power( dMultBase * dMultRawErr / (dMultRaw*dMultRaw),2.0);

    hScaled_Mult->SetBinContent(bin, dMultBase / dMultRaw);
    hScaled_Mult->SetBinError(bin, TMath::Sqrt(dErrSq));
  }

  TLegend* legMult = new TLegend(0.52,0.7,0.8,0.88);
  legMult->SetBorderSize(0);
  legMult->SetFillColorAlpha(0,0);
  // legMult->SetHeader(Form(""));
  legMult->AddEntry(hRaw_Mult,"p-Pb","p");
  legMult->AddEntry(hBase_MultInt,"pp","p");


  TCanvas* canMult = new TCanvas("canMult","canMult",800,400);
  canMult->Divide(2,1);
  canMult->cd(1);
  TH1* frame_Mult_1 = (TH1*) gPad->DrawFrame(hRaw_Mult->GetXaxis()->GetXmin(),0.0,hRaw_Mult->GetXaxis()->GetXmax(),1.2*hRaw_Mult->GetMaximum());
  frame_Mult_1->SetTitle("<M>; multiplicity %; <M>");
  hRaw_Mult->SetStats(0);
  hRaw_Mult->DrawCopy("same");
  hBase_MultInt->DrawCopy("same");
  legMult->Draw();
  canMult->cd(2);
  TH1* frame_Mult_2 = (TH1*) gPad->DrawFrame(hRaw_Mult->GetXaxis()->GetXmin(),0.0,hRaw_Mult->GetXaxis()->GetXmax(),1.3);
  frame_Mult_2->SetTitle("scaling factor k; multiplicity %; k");
  hScaled_Mult->DrawCopy("same");
  canMult->SaveAs(Form("%s/Mult.pdf",sOutPath.Data()));

  // TCanvas* canMult = new TCanvas("canMult","canMult",800,400);
  // canMult->Divide(2,1);
  // canMult->cd(1);
  // hRaw_Mult->SetStats(0);
  // hRaw_Mult->DrawCopy();
  // hBase_MultInt->DrawCopy("same");
  // canMult->cd(2);
  // hScaled_Mult->SetStats(0);
  // hScaled_Mult->DrawCopy();
  // canMult->SaveAs(Form("%s/Mult.pdf",sOutPath.Data()));


  // cn
  TH1D* hRaw_cn = LoadHisto(Form("hCum2_Refs_harm2_%s",sGapRaw.Data()),fileInRaw); if(!hRaw_cn) { return; }
  StyleHist(hRaw_cn, colRaw, markRaw);
  TH1D* hBase_cn_int = LoadHisto(Form("hCum2_Refs_harm2_%s",sGapBase.Data()),fileInBaseInt); if(!hBase_cn_int) { return; }
  StyleHist(hBase_cn_int, colBase, markBase);

  Double_t dBaseCn = hBase_cn_int->GetBinContent(1);
  Double_t dBaseCnErr = hBase_cn_int->GetBinError(1);

  TH1D* hSubt_cn = (TH1D*) hRaw_cn->Clone(Form("%s_subt",hRaw_cn->GetName()));  if(!hSubt_cn) { return; }
  StyleHist(hSubt_cn, colSubt, markSubt);

  TH1D* hBase_cn_scaled = (TH1D*) hRaw_cn->Clone(Form("%s_scaled",hRaw_cn->GetName()));  if(!hBase_cn_scaled) { return; }
  StyleHist(hBase_cn_scaled, colBase, markBaseScaled);

  for(Int_t cnBin(1); cnBin < hSubt_cn->GetNbinsX()+1; ++cnBin)
  {
    Double_t dFactor = hScaled_Mult->GetBinContent(cnBin);
    Double_t dFactorErr = hScaled_Mult->GetBinError(cnBin);

    Double_t dRawCn = hRaw_cn->GetBinContent(cnBin);
    Double_t dRawCnErr = hRaw_cn->GetBinError(cnBin);

    // subtraction
    Double_t dCont = dRawCn - dBaseCn * dFactor;
    Double_t dErrSq = TMath::Power(dRawCnErr,2.0) + TMath::Power(dBaseCnErr*dFactor,2.0) + TMath::Power(dBaseCn*dFactorErr,2.0);

    hSubt_cn->SetBinContent(cnBin, dCont);
    hSubt_cn->SetBinError(cnBin, TMath::Sqrt(dErrSq));

    // at this point, the subtracted cn is ready in hSubt_cn

    // scaling
    Double_t dContScaled = dBaseCn * dFactor;
    Double_t dErrScaled = TMath::Power(dBaseCnErr*dFactor,2.0) + TMath::Power(dBaseCn*dFactorErr,2.0);

    hBase_cn_scaled->SetBinContent(cnBin, dContScaled);
    hBase_cn_scaled->SetBinError(cnBin, TMath::Sqrt(dErrScaled));
  }

  TLegend* legRefs = new TLegend(0.12,0.6,0.4,0.88);
  legRefs->SetBorderSize(0);
  legRefs->SetFillColorAlpha(0,0);
  legRefs->SetHeader(Form("RFPs"));
  legRefs->AddEntry(hBase_cn_int,"pp (raw)","p");
  legRefs->AddEntry(hBase_cn_scaled,"pp (scaled)","p");
  legRefs->AddEntry(hRaw_cn,"p-Pb (raw)","p");
  legRefs->AddEntry(hSubt_cn,"p-Pb (subt)","p");


  TCanvas* canRef = new TCanvas("canRef","canRef", 600,600);
  canRef->cd(1);
  TH1* frame_Ref_1 = (TH1*) gPad->DrawFrame(hBase_cn_int->GetXaxis()->GetXmin(),-0.001+hSubt_cn->GetMinimum(),hBase_cn_int->GetXaxis()->GetXmax(),2.0*hBase_cn_int->GetMaximum());
  frame_Ref_1->SetTitle("c_{2}{2}; multiplicity %; c_{2}{2}");
  hBase_cn_int->DrawCopy("same");
  hRaw_cn->DrawCopy("same");
  hBase_cn_scaled->DrawCopy("same");
  hSubt_cn->DrawCopy("same");
  legRefs->Draw();
  canRef->SaveAs(Form("%s/Refs.pdf",sOutPath.Data()));

  // TCanvas* canRef = new TCanvas("canRef","canRef", 800,400);
  // canRef->Divide(2,1);
  // canRef->cd(1);
  // TH1* frame_Ref_1 = (TH1*) gPad->DrawFrame(hBase_cn_int->GetXaxis()->GetXmin(),-0.001+hSubt_cn->GetMinimum(),hBase_cn_int->GetXaxis()->GetXmax(),1.5*hBase_cn_int->GetMaximum());
  // frame_Ref_1->SetTitle("; multiplicity %; c_{2}{2}");
  // hBase_cn_int->DrawCopy("same");
  // hRaw_cn->DrawCopy("same");
  // hBase_cn_scaled->DrawCopy("same");
  // hSubt_cn->DrawCopy("same");
  // canRef->cd(2);
  // hSubt_cn->DrawCopy();
  // canRef->SaveAs(Form("%s/Refs.pdf",sOutPath.Data()));

  fileOut->cd();
  hRaw_Mult->Write();
  hBase_MultInt->Write();
  hRaw_cn->Write();
  hBase_cn_int->Write();
  hSubt_cn->Write();

  fileOutSubt->cd();
  hSubt_cn->Write();

  // dn
  for(Int_t iSpecies(0); iSpecies < iNumSpecies; ++iSpecies)
  {
    TString sSpecies = sSpecies_list[iSpecies];

    TH1D* hBase_dn_int = LoadHisto(Form("hCum2_%s_harm2_%s_cent0",sSpecies.Data(),sGapBase.Data()),fileInBaseInt); if(!hBase_dn_int) { return; }
    StyleHist(hBase_dn_int, colBase, markBase);

    TH1D* hBase_dn_int_scaled = (TH1D*) hBase_dn_int->Clone(Form("%s_scaled",hBase_dn_int->GetName()));
    StyleHist(hBase_dn_int_scaled, colBase, markBaseScaled);

    fileOut->cd();
    hBase_dn_int->Write(Form("hCum2_%s_harm2_%s_cent0_int",sSpecies.Data(),sGapBase.Data()));

    for(Int_t cnBin(1); cnBin < hSubt_cn->GetNbinsX()+1; ++cnBin)
    {
      Int_t cent = cnBin-1;
      TH1D* hRaw_dn = LoadHisto(Form("hCum2_%s_harm2_%s_cent%d",sSpecies.Data(),sGapRaw.Data(),cent),fileInRaw); if(!hRaw_dn) { return; }
      StyleHist(hRaw_dn, colRaw, markRaw);

      TH1D* hSubt_dn = (TH1D*) hRaw_dn->Clone(Form("%s_subt",hRaw_dn->GetName()));  if(!hSubt_dn) { return; }
      StyleHist(hSubt_dn, colSubt, markSubt);

      Double_t dFactor = hScaled_Mult->GetBinContent(cnBin);
      Double_t dFactorErr = hScaled_Mult->GetBinError(cnBin);

      for(Int_t pt(1); pt < hSubt_dn->GetNbinsX()+1; ++pt)
      {
        Double_t dRawDn = hRaw_dn->GetBinContent(pt);
        Double_t dRawDnErr = hRaw_dn->GetBinError(pt);

        Double_t dBaseDn = hBase_dn_int->GetBinContent(pt);
        Double_t dBaseDnErr = hBase_dn_int->GetBinError(pt);

        // subtraction

        Double_t dCont = dRawDn - dBaseDn * dFactor;
        Double_t dErrSq = TMath::Power(dRawDnErr,2.0) + TMath::Power(dBaseDnErr*dFactor,2.0) + TMath::Power(dBaseDn*dFactorErr,2.0);

        hSubt_dn->SetBinContent(pt,dCont);
        hSubt_dn->SetBinError(pt,TMath::Sqrt(dErrSq));

        // scaling
        Double_t dContScaled = dBaseDn * dFactor;
        Double_t dErrScaled = TMath::Power(dBaseDnErr*dFactor,2.0) + TMath::Power(dBaseDn*dFactorErr,2.0);
        hBase_dn_int_scaled->SetBinContent(pt,dContScaled);
        hBase_dn_int_scaled->SetBinError(pt,TMath::Sqrt(dErrScaled));
      }

      // making vn^subt out of cn^sub and dn^sub
      TH1D* hSubt_vn = (TH1D*) hSubt_dn->Clone(Form("hFlow2_%s_harm2_%s_cent%d",sSpecies.Data(), sGap.Data(), cent));
      hSubt_vn->Scale(1.0/TMath::Sqrt(hSubt_cn->GetBinContent(cnBin)));

      TLegend* legDiff = new TLegend(0.12,0.6,0.4,0.88);
      legDiff->SetBorderSize(0);
      legDiff->SetFillColorAlpha(0,0);
      legDiff->SetHeader(Form("%s (%s)",sSpecies.Data(), sCentLabel[cent].Data()));
      legDiff->AddEntry(hBase_dn_int,"pp (raw)","p");
      legDiff->AddEntry(hBase_dn_int_scaled,"pp (scaled)","p");
      legDiff->AddEntry(hRaw_dn,"p-Pb (raw)","p");
      legDiff->AddEntry(hSubt_dn,"p-Pb (subt)","p");

      TCanvas* canDiff = new TCanvas("canDiff","canDiff",600,600);
      canDiff->cd(1);
      TH1* frame_Diff1 = (TH1*) gPad->DrawFrame(hBase_dn_int->GetXaxis()->GetXmin(),-0.001+hSubt_dn->GetMinimum(),hBase_dn_int->GetXaxis()->GetXmax(),1.5*hBase_dn_int->GetMaximum());
      frame_Diff1->SetTitle("d_{2}{2}; p_{T} (GeV/c); d_{2}{2}");
      hBase_dn_int->SetStats(0);
      hBase_dn_int->SetMinimum(0.0);
      hBase_dn_int->DrawCopy("same");
      hBase_dn_int_scaled->DrawCopy("same");
      hRaw_dn->DrawCopy("same");
      hSubt_dn->DrawCopy("same");
      legDiff->Draw();
      canDiff->SaveAs(Form("%s/%s_cent%d.pdf",sOutPath.Data(),sSpecies.Data(),cent));

      // TCanvas* canDiff = new TCanvas("canDiff","canDiff",800,400);
      // canDiff->Divide(2,1);
      // canDiff->cd(1);
      // hBase_dn_int->SetMinimum(0.0);
      // hBase_dn_int->DrawCopy();
      // hBase_dn_int_scaled->DrawCopy("same");
      // hRaw_dn->DrawCopy("same");
      // hSubt_dn->DrawCopy("same");
      // canDiff->cd(2);
      // hSubt_dn->DrawCopy();
      // canDiff->SaveAs(Form("%s/%s_cent%d.pdf",sOutPath.Data(),sSpecies.Data(),cent));

      // save
      fileOutSubt->cd();
      hSubt_vn->Write();

      fileOut->cd();
      hBase_dn_int->Write();
      hRaw_dn->Write();
      hSubt_dn->Write();
    }
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
TH1D* Scale(TH1D* base, Double_t factor, Double_t factorErrSq)
{
  if(!base) { printf("ERROR-Scale: Hist 'base' does not found.\n"); return 0x0; }

  TH1D* scale = (TH1D*) base->Clone(Form("%s_scaled",base->GetName()));
  if(!scale) { printf("ERROR-scale: Hist 'scale' does not cloned properly.\n"); return 0x0; }

  for(Int_t bin(1); bin < base->GetNbinsX()+1; ++bin)
  {
    Double_t dCon = base->GetBinContent(bin);
    Double_t dErr = base->GetBinError(bin);

    Double_t dErrSq = factorErrSq*dCon*dCon + factor*factor*dErr*dErr;

    scale->SetBinContent(bin, factor * dCon);
    scale->SetBinError(bin, TMath::Sqrt(dErrSq));
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
TH1D* Subtract(TH1D* raw, TH1D* base)
{
  if(!raw) { printf("ERROR-Subtract: Hist 'raw' does not found.\n"); return 0x0; }
  if(!base) { printf("ERROR-Subtract: Hist 'base' does not found.\n"); return 0x0; }

  TH1D* sub = (TH1D*) raw->Clone(Form("%s_sub",raw->GetName()));
  if(!sub) { printf("ERROR-Subtract: Hist 'sub' does not cloned properly.\n"); return 0x0; }

  for(Int_t bin(1); bin < sub->GetNbinsX()+1; ++bin)
  {
    Double_t con_raw = raw->GetBinContent(bin);
    Double_t err_raw = raw->GetBinError(bin);
    Double_t con_base = base->GetBinContent(bin);
    Double_t err_base = base->GetBinError(bin);

    // sub->SetBinContent(bin, con_raw - factor * con_base);
    // sub->GetBinError(bin, 0.01*con_raw);
    // sub->GetBinError(bin, err_raw*err_raw + factor*factor*err_base*err_base );

  }
  return sub;
}
// ==================================================================================================================
TH1D* Subtract_new_cn(TH1D* raw, TH1D* base, TProfile* mult_raw, TProfile* mult_base)
{
  if(!raw) { printf("ERROR-Subtract_new_cn: Hist 'raw' does not found.\n"); return 0x0; }
  if(!base) { printf("ERROR-Subtract_new_cn: Hist 'base' does not found.\n"); return 0x0; }
  if(!mult_raw) { printf("ERROR-Subtract_new_cn: Hist 'mult_raw' does not found.\n"); return 0x0; }
  if(!mult_base) { printf("ERROR-Subtract_new_cn: Hist 'mult_base' does not found.\n"); return 0x0; }

  TH1D* sub = (TH1D*) raw->Clone(Form("%s_sub",raw->GetName()));
  if(!sub) { printf("ERROR-Subtract_new_cn: Hist 'sub' does not cloned properly.\n"); return 0x0; }

  Double_t dMult_base = mult_base->GetBinContent(1);
  Double_t dMult_base_err = mult_base->GetBinError(1);
  Double_t con_base = base->GetBinContent(1);
  Double_t err_base = base->GetBinError(1);

  for(Int_t bin(1); bin < sub->GetNbinsX()+1; ++bin)
  {
    Double_t dMult_raw = mult_raw->GetBinContent(bin);
    Double_t dMult_raw_err = mult_raw->GetBinError(bin);

    Double_t con_raw = raw->GetBinContent(bin);
    Double_t err_raw = raw->GetBinError(bin);

    if(dMult_base <= 0.0) { printf("ERROR-Subtract_new_cn dMultBase is 0.\n"); return 0x0; }
    Double_t dFactor = dMult_base / dMult_raw;
    Double_t dFactor_err_sq = TMath::Power(dMult_base_err/dMult_raw,2.0) + TMath::Power(dMult_base*dMult_raw_err,2.0) * TMath::Power(dMult_raw,-4.0);

    Double_t dCon = con_raw - dFactor*con_base;
    Double_t dErr = err_raw*err_raw + err_base*err_base*dFactor*dFactor + con_base*con_base*dFactor_err_sq;

    sub->SetBinContent(bin, dCon);
    sub->GetBinError(bin, TMath::Sqrt(dErr));

  }
  return sub;
}
// ==================================================================================================================
TH1D* Subtract_new_dn(TH1D* raw, TH1D* base, Double_t dFactor, Double_t dFactor_err_sq)
{
  if(!raw) { printf("ERROR-Subtract_new_dn: Hist 'raw' does not found.\n"); return 0x0; }
  if(!base) { printf("ERROR-Subtract_new_dn: Hist 'base' does not found.\n"); return 0x0; }
  // if(!mult_raw) { printf("ERROR-Subtract_new_dn: Hist 'mult_raw' does not found.\n"); return 0x0; }
  // if(!mult_base) { printf("ERROR-Subtract_new_dn: Hist 'mult_base' does not found.\n"); return 0x0; }

  TH1D* sub = (TH1D*) raw->Clone(Form("%s_sub",raw->GetName()));
  if(!sub) { printf("ERROR-Subtract_new_dn: Hist 'sub' does not cloned properly.\n"); return 0x0; }

  // Double_t dMult_base = mult_base->GetBinContent(1);
  // Double_t dMult_base_err = mult_base->GetBinError(1);

  for(Int_t bin(1); bin < sub->GetNbinsX()+1; ++bin)
  {
    // Double_t dMult_raw = mult_raw->GetBinContent(bin);
    // Double_t dMult_raw_err = mult_raw->GetBinError(bin);

    Double_t con_raw = raw->GetBinContent(bin);
    Double_t err_raw = raw->GetBinError(bin);
    Double_t con_base = base->GetBinContent(bin);
    Double_t err_base = base->GetBinError(bin);

    // if(dMult_raw <= 0.0) { printf("ERROR-Subtract_new_dn: dMultBase is 0.\n"); return 0x0; }
    // Double_t dFactor = dMult_base / dMult_raw;
    // Double_t dFactor_err_sq = TMath::Power(dMult_base_err/dMult_raw,2.0) + TMath::Power(dMult_base*dMult_raw_err,2.0) * TMath::Power(dMult_raw,-4.0);

    Double_t dCon = con_raw - dFactor*con_base;
    Double_t dErr = err_raw*err_raw + err_base*err_base*dFactor*dFactor + con_base*con_base*dFactor_err_sq;

    sub->SetBinContent(bin, dCon);
    sub->GetBinError(bin, TMath::Sqrt(dErr));

  }
  return sub;
}
// ==================================================================================================================
