#include "TFile.h"
#include "TH1.h"
#include "TList.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"
#include "TSystem.h"

TH1D* DivideHistos(TH1* nom, TH1* denom, Bool_t bCor = kFALSE);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle, Bool_t showStats = kFALSE);
TList* PrepareCanvas(TCanvas* can);


void CompareRun1()
{
  TString sInputFileHEP = "/Users/vpacik/NBI/Flow/data/HEPdata/PLB726(2013)/extracted-v2SP.root";
  TString sInputFileNew = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/merged-pPb-16qt-nua/output_run1comp/gap08/Processed.root";
  TString sOutputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/outcome/compRun1/merged-pPb-16qt-nua/output_run1comp/";

  TString sCentLabels[] = {"0-20%","20-40%","40-60%","60-100%"};
  Int_t iNumCent = sizeof(sCentLabels)/sizeof(sCentLabels[0]);

  Int_t iNumSpecies = 4;
  TString sSpecies[] = {"Charged","Pion","Kaon","Proton"};
  Color_t colors[] = {kBlack,kRed,kBlue,kGreen+2};
  Int_t symbolsHEP[] = {kOpenSquare,kOpenSquare,kOpenSquare,kOpenSquare};
  Int_t symbolsNew[] = {kFullCircle,kFullCircle,kFullCircle,kFullCircle};


  TFile* fileInputHEP = TFile::Open(sInputFileHEP.Data(),"READ");
  if(!fileInputHEP) { printf("ERROR : fileInputHEP not found!\n"); return; }
  TFile* fileInputNew = TFile::Open(sInputFileNew.Data(),"READ");
  if(!fileInputNew) { printf("ERROR : fileInputNew not found!\n"); return; }

  fileInputHEP->ls();
  fileInputNew->ls();



  for(Int_t iCent(0); iCent < iNumCent; ++iCent)
  {
    TCanvas* canCent = new TCanvas("canCent","canCent",500,1200);
    TList* listCent = PrepareCanvas(canCent);
    if(!listCent) { printf("ERROR : PrepareCanvas (CEtn) failed\n"); return; }
    TPad* padMainCent = (TPad*) listCent->At(0);
    TPad* padRatioCent = (TPad*) listCent->At(1);
    TH1* frameMainCent = (TH1*) listCent->At(2);
    TH1* frameRatioCent = (TH1*) listCent->At(3);
    frameRatioCent->SetTitle(Form("; p_{T} (GeV/c); Run2/Run1   "));

    TLegend* legCent = new TLegend(0.2,0.7,0.4,0.89);
    legCent->SetFillColorAlpha(0,0);
    legCent->SetBorderSize(0);
    legCent->SetHeader(Form("(%s)",sCentLabels[iCent].Data()));


    for(Int_t iSpecies(0); iSpecies < iNumSpecies; ++iSpecies)
    {
      TString sSpeciesName = sSpecies[iSpecies];
      TH1F* hHEP = (TH1F*) fileInputHEP->Get(Form("hFlow2_%s_harm2_gap08_cent%d",sSpeciesName.Data(),iCent));
      if(!hHEP) { printf("ERROR : hHEP 'hFlow2_%s_harm2_gap08_cent%d' not found!\n",sSpeciesName.Data(),iCent); return; }

      TH1D* hNew = (TH1D*) fileInputNew->Get(Form("hFlow2_%s_harm2_gap08_cent%d",sSpeciesName.Data(),iCent));
      if(!hNew) { printf("ERROR : hNew 'hFlow2_%s_harm2_gap08_cent%d' not found!\n",sSpeciesName.Data(),iCent); return; }

      StyleHist(hNew,colors[iSpecies],symbolsNew[iSpecies]);
      StyleHist(hHEP,colors[iSpecies],symbolsHEP[iSpecies]);

      TH1D* hRatio = DivideHistos(hNew,hHEP,kFALSE);

      TLegend* leg = new TLegend(0.2,0.7,0.4,0.89);
      leg->SetFillColorAlpha(0,0);
      leg->SetBorderSize(0);
      leg->SetHeader(Form("%s (%s)",sSpeciesName.Data(),sCentLabels[iCent].Data()));
      leg->AddEntry(hNew,"v_{2}{2} (Run2)","pl");
      leg->AddEntry(hHEP,"v_{2}{2,SP} (Run1)","pl");

      TCanvas* can = new TCanvas(Form("can"),Form("can"),500,1200);
      TList* list = PrepareCanvas(can);
      if(!list) { printf("ERROR : PrepareCanvas failed\n"); return; }
      // list->ls();
      TPad* padMain = (TPad*) list->At(0);
      TPad* padRatio = (TPad*) list->At(1);
      TH1* frameMain = (TH1*) list->At(2);
      TH1* frameRatio = (TH1*) list->At(3);

      frameRatio->SetTitle(Form("; p_{T} (GeV/c); Run2/Run1   "));

      padMain->cd();
      hNew->DrawCopy("same");
      hHEP->DrawCopy("same");
      leg->Draw();

      padRatio->cd();
      hRatio->Draw("same");

      padMainCent->cd();
      hNew->DrawCopy("same");
      hHEP->DrawCopy("same");

      padRatioCent->cd();
      hRatio->Draw("same");

      legCent->AddEntry(hNew,Form("v_{2}{2} %s (Run2)",sSpeciesName.Data()),"pl");
      legCent->AddEntry(hHEP,Form("v_{2}{2,SP} %s (Run1)",sSpeciesName.Data()),"pl");

      gSystem->mkdir(sOutputPath.Data(),kTRUE);
      can->SaveAs(Form("%s/compRun1_%s_cent%d.pdf",sOutputPath.Data(),sSpeciesName.Data(),iCent),"pdf");
    }
    padMainCent->cd();
    legCent->Draw();
    canCent->SaveAs(Form("%s/compRun1_ALL_cent%d.pdf",sOutputPath.Data(),iCent),"pdf");
  }

  return;
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
TH1D* DivideHistos(TH1* nom, TH1* denom, Bool_t bCor)
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
  for(Short_t iBinNom(1); iBinNom < binsNom+1; iBinNom++)
  {
    Double_t dBinCenter = nom->GetBinCenter(iBinNom);
    Int_t iBinDenom = denom->FindFixBin(dBinCenter);

    dContNom = nom->GetBinContent(iBinNom);
    dErrNom = nom->GetBinError(iBinNom);
    dContDenom = denom->GetBinContent(iBinDenom);
    dErrDenom = denom->GetBinError(iBinDenom);

    if(dContDenom == 0.0) continue;

    dContRatio =  dContNom / dContDenom;
    dErrRatio = TMath::Power(dErrNom/dContDenom, 2) + TMath::Power( dErrDenom*dContNom/(dContDenom*dContDenom), 2);
    // printf("Err (before) : %g | ", TMath::Sqrt(dErrRatio));

    if(bCor) dErrRatio -= (2*dContNom*dErrDenom*dErrNom/TMath::Power(dContDenom,3));
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    ratio->SetBinContent(iBinNom,dContRatio);
    ratio->SetBinError(iBinNom,TMath::Sqrt(dErrRatio));
  }

  return ratio;
}
// ==================================================================================================================
TList* PrepareCanvas(TCanvas* can)
{
  if(!can) { printf("ERROR-PrepareCanvas: no canvas found!\n"); return 0x0; }

  TLine* lUnity = new TLine();
  lUnity->SetLineColor(kGray+1);
  lUnity->SetLineStyle(kDashed);

  can->cd();
  TPad* padMain = new TPad("padMain","padMain", 0, 0.3, 1, 1.0);
  padMain->SetBottomMargin(0.0);
  padMain->SetRightMargin(0.03);
  padMain->SetLeftMargin(0.13);
  padMain->Draw();
  padMain->cd();
  TH1* frame_canDiff_1 = (TH1*) gPad->DrawFrame(0.0,-0.05,4.0,0.4,Form("; p_{T} (GeV/c); v_{2}{2,|#Delta#eta|>0.8}"));
  frame_canDiff_1->SetTitleFont(43,"XY");
  frame_canDiff_1->SetTitleSize(18,"XY");
  frame_canDiff_1->SetTitleOffset(4.3,"X");
  frame_canDiff_1->SetLabelFont(43,"X");
  frame_canDiff_1->SetLabelSize(18,"X");
  // frame_canDiff_1->SetTitleFont(43,"Y");
  frame_canDiff_1->SetTitleOffset(2.2,"Y");

  can->cd();
  TPad* padRatio = new TPad("padRatio","padRatio", 0, 0.0, 1, 0.3);
  padRatio->SetTopMargin(0.0);
  padRatio->SetBottomMargin(0.25);
  padRatio->SetRightMargin(0.03);
  padRatio->SetLeftMargin(0.13);
  padRatio->Draw();
  padRatio->cd();
  TH1* frame_canDiff_2 = (TH1*) gPad->DrawFrame(0.0,0.65,4.0,1.35);
  frame_canDiff_2->SetTitle(Form("; p_{T} (GeV/c); Kch / K0s   "));
  frame_canDiff_2->SetNdivisions(505,"Y");
  frame_canDiff_2->SetTitleFont(43,"XY");
  frame_canDiff_2->SetTitleSize(18,"XY");
  frame_canDiff_2->SetTitleOffset(4.3,"X");
  frame_canDiff_2->SetTitleOffset(2.2,"Y");
  frame_canDiff_2->SetLabelFont(43,"XY");
  frame_canDiff_2->SetLabelSize(18,"XY");
  lUnity->DrawLine(0.0,1.0,4.0,1.0);

  TList* list = new TList();
  list->Add(padMain);
  list->Add(padRatio);
  list->Add(frame_canDiff_1);
  list->Add(frame_canDiff_2);

  return list;
}
