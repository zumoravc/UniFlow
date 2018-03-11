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

TH1D* Scale(TH1D* base, Double_t factor);
TH1D* ScaleCn(TH1D* base, TH1D* mult);
TH1D* Subtract(TH1D* raw, TH1D* base, Double_t factor = 1.0);

// colors for centrality
Color_t colors[] = {kGreen+2, kBlue, kBlack, kMagenta+1};


void Subt_ppb_pp(TString sGap="gap00")
{
  TString sSpecies = "Proton";
  // TString sMethod = "GF_eventweighted";
  // TString sOutputTag = "output_vn";
  // TString sOutputTagInt = sOutputTag + "_int";

  TString sGapBase = sGap;
  TString sGapRaw = sGapBase;

  TString sInFileRaw = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pPb-16qt/output/" + sGapBase; // + "/" + sMethod;
  TString sInFileBaseInt = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pp-16kl/output_int/" + sGapBase; // + "/"; + sMethod;
  // TString sOutFolder = sInFileRaw+"/"+sMethod+"/pPb_pp_subt_"+sGapRaw+"/"+sSpecies;
  TString sOutFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/pPb-16qt/subt-pPb-pp/"+sGapBase+"/"+sSpecies;
  TString sOutFile = sOutFolder+"/Subt_results.root";

  const Int_t iNumCent = 4;
  TString sCentLabel[iNumCent] = {"0-20%", "20-40%", "40-60%", "60-100%"};

  // ==================================================================================================================

  // === LOADING INPUT ===
  // multiplicities
  TFile* fileInRaw_Mult = OpenFile(sInFileRaw+"/Mult.root"); if(!fileInRaw_Mult) { return; }
  TProfile* hRaw_Mult = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInRaw_Mult); if(!hRaw_Mult) { return; }

  TFile* fileInBase_MultInt = OpenFile(sInFileBaseInt+"/Mult.root"); if(!fileInBase_MultInt) { return; }
  TProfile* hBase_MultInt = (TProfile*) LoadHisto("fpRefsMult_rebin",fileInBase_MultInt); if(!hBase_MultInt) { return; }
  Double_t dMult_Base_Int = hBase_MultInt->GetBinContent(1);

  // output files
  gSystem->mkdir(sOutFolder,kTRUE);  // Making output folder
  TFile* fileOut = OpenFile(sOutFile,"RECREATE"); if(!fileOut) { return; }

  // input files
  TFile* fileInRaw = OpenFile(sInFileRaw+"/Processed.root"); if(!fileInRaw) { return; }
  TFile* fileInBaseInt = OpenFile(sInFileBaseInt+"/Processed.root"); if(!fileInBaseInt) { return; }

  // cn
  TH1D* hRaw_cn = LoadHisto(Form("hCum2_Refs_harm2_%s",sGapRaw.Data()),fileInRaw); if(!hRaw_cn) { return; }
  StyleHist(hRaw_cn, kRed, kOpenSquare);
  TH1D* hBase_cn_int = LoadHisto(Form("hCum2_Refs_harm2_%s",sGapBase.Data()),fileInBaseInt); if(!hBase_cn_int) { return; }
  StyleHist(hBase_cn_int, kBlue, kFullCircle);

  TH1D* hRaw_cn_scaled = ScaleCn(hRaw_cn,hRaw_Mult);
  TH1D* hBase_cn_int_scaled = ScaleCn(hBase_cn_int,hBase_MultInt);

  // dn
  TList* list_Raw_dn = new TList();
  TList* list_Raw_dn_scaled = new TList();

  for(Int_t cent(0); cent < iNumCent; ++cent)
  {
    TH1D* temp = LoadHisto(Form("hCum2_%s_harm2_%s_cent%d",sSpecies.Data(),sGapRaw.Data(),cent),fileInRaw); if(!temp) { return; }
    StyleHist(temp, kRed, kOpenSquare);
    list_Raw_dn->Add(temp);

    TH1D* raw_scaled = Scale((TH1D*) list_Raw_dn->At(cent), hRaw_Mult->GetBinContent(cent+1)); if(!raw_scaled) { return; }
    list_Raw_dn_scaled->Add(raw_scaled);
  }

  TH1D* hBase_dn_int = LoadHisto(Form("hCum2_%s_harm2_%s_cent0",sSpecies.Data(),sGapBase.Data()),fileInBaseInt); if(!hBase_dn_int) { return; }
  StyleHist(hBase_dn_int, kBlue, kFullCircle);
  TH1D* hBase_dn_int_scaled = Scale(hBase_dn_int, hBase_MultInt->GetBinContent(1));

  fileOut->cd();
  hRaw_cn->Write("hRaw_cn");
  hBase_cn_int->Write("hBase_cn_int");
  hRaw_cn_scaled->Write("hRaw_cn_scaled");
  hBase_cn_int_scaled->Write("hBase_cn_int_scaled");

  list_Raw_dn->Write("list_Raw_dn",TObject::kSingleKey);
  hBase_dn_int->Write("hBase_dn_int");
  list_Raw_dn_scaled->Write("list_Raw_dn_scaled",TObject::kSingleKey);
  hBase_dn_int_scaled->Write("hBase_dn_int_scaled");

  // === SUBTRACTING pPb(cent) - pp (MB) ===

  // cn{2}^sub = <M>^raw * cn{2}^raw - <M>^base * cn{2}^base
  TH1D* hSub_cn_int = Subtract(hRaw_cn_scaled, hBase_cn_int_scaled,1.0); if(!hSub_cn_int) { return; }

  // dn{2}^sub = <M>dn{2}^raw - <M>dn{2}^base
  TList* list_SubPP_dn_int = new TList();
  for(Int_t centRaw(0); centRaw < iNumCent; ++centRaw)
  {
    TH1D* temp_raw_cent = (TH1D*) list_Raw_dn_scaled->At(centRaw);

    // subtracting (int) pp
    TH1D* hSubPP_dn_int = Subtract(temp_raw_cent, hBase_dn_int_scaled);
    list_SubPP_dn_int->Add(hSubPP_dn_int);
  }

  // Making vn^sub out of dn^sub / sqrt(cn^sub)
  TList* list_SubPP_vn_int = new TList();
  for(Int_t centRaw(0); centRaw < iNumCent; ++centRaw)
  {
    TH1D* hSub_dn = (TH1D*) list_SubPP_dn_int->At(centRaw); if(!hSub_dn) { return; }

    // dividing by cn{2}^sub (int)
    TH1D* hSubPP_vn_int = (TH1D*) hSub_dn->Clone(Form("%s_vn_int", hSub_dn->GetName())); if(!hSubPP_vn_int) { return; }
    StyleHist(hSubPP_vn_int, kBlue, kFullCircle);
    hSubPP_vn_int->Scale(1.0/TMath::Sqrt(hSub_cn_int->GetBinContent(1)));
    hSubPP_vn_int->SetName(Form("hSubPP_vn_int_cent%d",centRaw));
    list_SubPP_vn_int->Add(hSubPP_vn_int);

  }

  fileOut->cd();
  hSub_cn_int->Write("hSub_cn_int");
  list_SubPP_dn_int->Write("list_SubtPP_dn_ppint",TObject::kSingleKey);
  list_SubPP_vn_int->Write("list_SubtPP_vn_ppint",TObject::kSingleKey);

  // UNIVERSALL PLOTTING
  // RFPs
  TLegend* leg_Refs = new TLegend(0.12,0.12,0.6,0.3);
  leg_Refs->SetBorderSize(0.);
  leg_Refs->SetFillColor(0);
  leg_Refs->AddEntry(hRaw_cn,"pPb","p");
  leg_Refs->AddEntry(hBase_cn_int,"pp (0-100%)","p");

  TCanvas* canRefs = new TCanvas("canRefs","canRefs",1200,400);
  canRefs->Divide(3,1);
  canRefs->cd(1);
  TH1* frame_Ref = (TH1*) gPad->DrawFrame(0,0,100,0.01);
  frame_Ref->SetTitle("raw <<2>>; cent %");
  hRaw_cn->Draw("same");
  hBase_cn_int->Draw("same");
  leg_Refs->Draw();
  canRefs->cd(2);
  TH1* frame_Ref_2 = (TH1*) gPad->DrawFrame(0,-1.0,100,10.0);
  frame_Ref_2->SetTitle("<M>^{2} * <<2>>; cent %");
  hRaw_cn_scaled->Draw("same");
  hBase_cn_int_scaled->Draw("same");
  canRefs->cd(3);
  TH1* frame_Ref_3 = (TH1*) gPad->DrawFrame(0,0,100,10.0);
  frame_Ref_3->SetTitle("<M>^{raw,2} <<2>> - <M>^{base,2}<<2>>; cent %");
  hSub_cn_int->Draw("same");
  canRefs->SaveAs(Form("%s/cn_subt.pdf",sOutFolder.Data()),"pdf");

  for(Int_t centRaw(0); centRaw < iNumCent; ++centRaw)
  {
    TH1D* temp_raw_cent = (TH1D*) list_Raw_dn_scaled->At(centRaw);
    TH1D* hSubPP_dn_int = (TH1D*) list_SubPP_dn_int->At(centRaw);


    TLegend* leg = new TLegend(0.12,0.5,0.6,0.89);
    leg->SetBorderSize(0.);
    leg->SetFillColor(0);
    leg->AddEntry(temp_raw_cent,Form("pPb (unsub, %s)",sCentLabel[centRaw].Data()),"p");
    leg->AddEntry(hBase_dn_int,"pp (0-100%)","p");

    TCanvas* can = new TCanvas("can","can",1200,400);
    can->Divide(3,1);
    can->cd(1);
    TH1* frame = gPad->DrawFrame(0,0,10,0.2);
    frame->SetTitle("raw <<2'>>; p_{T} (GeV/c)");
    ((TH1D*) list_Raw_dn->At(centRaw))->Draw("same");
    hBase_dn_int->Draw("same");

    leg->Draw();

    can->cd(2);
    TH1* frame2 = gPad->DrawFrame(0,-0.03,10,1.0);
    frame2->SetTitle("<M> * <<2'>>; p_{T} (GeV/c)");
    temp_raw_cent->Draw("same");
    hBase_dn_int_scaled->Draw("same");

    can->cd(3);
    TH1* frame3 = gPad->DrawFrame(0,-0.05,10,1.0);
    frame3->SetTitle("<M>^{pPb}<<2'>>^{pPb} - <M>^{pp}<<2'>>^{pp}; p_{T} (GeV/c)");
    hSubPP_dn_int->Draw("same");


    can->SaveAs(Form("%s/Subt_pp-pbp_cent%d.pdf",sOutFolder.Data(),centRaw),"pdf");
  }
  for(Int_t centRaw(0); centRaw < iNumCent; ++centRaw)
  {
    TH1D* hSubPP_vn_int = (TH1D*) list_SubPP_vn_int->At(centRaw);

    TLegend* leg = new TLegend(0.12,0.65,0.6,0.89);
    leg->SetBorderSize(0.);
    leg->SetFillColor(0);
    leg->AddEntry(hSubPP_vn_int,"pp (0-100%)","p");

    TCanvas* canVn = new TCanvas("canVn","canVn",400,400);
    // canVn->Divide(,1);
    canVn->cd(1);
    TH1* frame_vn = (TH1*) gPad->DrawFrame(0.,-0.05,10.,0.3);
    frame_vn->SetTitle(Form("v_{2}{2}^{sub} (%s V0A); p_{T} (GeV/c)",sCentLabel[centRaw].Data()));
    hSubPP_vn_int->Draw("same");

    // canVn->cd(2);
    // canVn->cd(3);
    leg->Draw();

    canVn->SaveAs(Form("%s/Subt_ppb_pp_vn_cent%d.pdf",sOutFolder.Data(),centRaw),"pdf");
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
TH1D* Subtract(TH1D* raw, TH1D* base, Double_t factor)
{
  if(!raw) { printf("ERROR-Subtract: Hist 'raw' does not found.\n"); return 0x0; }
  if(!base) { printf("ERROR-Subtract: Hist 'base' does not found.\n"); return 0x0; }

  TH1D* sub = (TH1D*) base->Clone(Form("%s_sub",raw->GetName()));
  if(!sub) { printf("ERROR-Subtract: Hist 'sub' does not cloned properly.\n"); return 0x0; }

  for(Int_t bin(1); bin < sub->GetNbinsX()+1; ++bin)
  {
    Double_t con_raw = raw->GetBinContent(bin);
    Double_t err_raw = raw->GetBinError(bin);
    Double_t con_base = base->GetBinContent(bin);
    Double_t err_base = base->GetBinError(bin);

    sub->SetBinContent(bin, con_raw - factor * con_base);
    // sub->GetBinError(bin, 0.01*con_raw);
    sub->GetBinError(bin, err_raw*err_raw + factor*factor*err_base*err_base );

  }
  return sub;
}
// ==================================================================================================================
