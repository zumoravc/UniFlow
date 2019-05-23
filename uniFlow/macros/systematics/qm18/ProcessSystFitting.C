#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TList.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TMath.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TColor.h"

TH1D* DivideHistos(TH1* nom, TH1* denom, Bool_t bCor = kFALSE);
TH1D* BarlowTest(TH1* nom, TH1* denom, Bool_t bCor = kFALSE);
TH1D* ApplyBarlow(TH1* diff, TH1* barlow, Double_t dCut);
TH1D* Diff(TH1* ratio);
void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle, Bool_t showStats = kFALSE);
TList* PrepareCanvas(TCanvas* can);
void SetCustomPalette();

// void ProcessSyst(Int_t iCent, Int_t iTag, Int_t iSpecies, Bool_t bCorrelated);
TH1D* GetRMS(TList* list);

Color_t colBase = kRed;
Color_t colSyst = kBlue+2;

Int_t markBase = kOpenSquare;
Int_t markSyst = kFullCircle;

Double_t dRatioYmin = 0.7;
Double_t dRatioYmax = 1.3;

Double_t dDiffYmin = 0.0;
Double_t dDiffYmax = 0.2;

Double_t dBarlowCut = -1.0;

Double_t dFitXmin = 0.2;
Double_t dFitXmax = 7.0;
// Double_t dFitXmax = 1.2;

Int_t iCent = 3;
Int_t iTag = 0;
Int_t iSpecies = 0;

Bool_t bCorrelated = 0;
Bool_t bCorRatio = 0;

TString sFitFunc = "[0]"; Int_t iNumParFit = 1;
// TString sFitFunc = "[0]+[1]*x+[2]*x*x"; Int_t iNumParFit = 3;

TString sGap = "gap04";

TString sCentLabels[] = {"0-10%","10-20%","20-40%","40-60%"};

TString sType = "fitting"; TString sTags[] = {"flowBG","massSig","massBG","range1","range2","range3","range4","range5","range6"}; TString sSpeciesList[] = {"K0s"};
// TString sType = "fitting"; TString sTags[] = {"flowBG","massSig","massBG","range1","range2","range3","range4"}; TString sSpeciesList[] = {"Lambda"};
// TString sType = "fitting"; TString sTags[] = {"massBG","massSig","range1","range2","range3","range4"}; TString sSpeciesList[] = {"Phi"};

TString sInputSyst = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/out-binning-3/"+sType;
TString sOutput = sInputSyst;

TFile* fileOutput = 0x0;

Int_t iNumTags = sizeof(sTags)/sizeof(sTags[0]);
Int_t iNumCent = sizeof(sCentLabels)/sizeof(sCentLabels[0]);
Int_t iNumSpecies = sizeof(sSpeciesList)/sizeof(sSpeciesList[0]);

void ProcessSystFitting()
{
  TString sSpecies = sSpeciesList[0];

  SetCustomPalette();

  fileOutput = TFile::Open(Form("%s/fitting_%s.root",sOutput.Data(),sSpecies.Data()),"RECREATE");
  if(!fileOutput) { printf("no fileOutput\n"); return; }

  for(Int_t iCent(0); iCent< iNumCent; ++iCent)
  {
    TList* listSyst = new TList();

    TLegend* leg = new TLegend(0.45,0.3,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0,0);

    for(Int_t tag(0); tag < iNumTags; ++tag)
    {
      TString sTag = sTags[tag];
      TFile* fileInput = TFile::Open(Form("%s/%s/%s.root",sInputSyst.Data(), sSpecies.Data(), sTag.Data()),"READ");
      if(!fileInput)  { printf("ERROR : fileInput not found!\n"); return; }

      // fileInput->ls();
      TList* list = (TList*) fileInput->Get(Form("%s_cent%d",sSpecies.Data(),iCent));
      // list->ls();
      TH1D* hist = (TH1D*) list->FindObject("histRatio")->Clone(Form("histRatio_%s",sTag.Data()));
      if(!hist) { printf("ERROR : hist not found!\n"); return; }

      listSyst->Add(hist);
      leg->AddEntry(hist,sTag.Data(),"pl");
      // fileInput->Close();
    }


    listSyst->ls();
    printf("Entries %d \n", listSyst->GetEntries());
    TH1D* histSTD = GetRMS(listSyst);
    leg->AddEntry(histSTD,"Mean (RMS)","pl");

    TLatex* latex = new TLatex();
    latex->SetTextSize(0.05);
    latex->SetNDC();

    TCanvas* can = new TCanvas("can","can",800,400);
    can->Divide(2,1);
    can->cd(1);
    TH1* frame = gPad->DrawFrame(histSTD->GetXaxis()->GetXmin(), 0.85,histSTD->GetXaxis()->GetXmax(), 1.15);
    frame->SetTitle("x_{syst}/x_{default}; p_{T} (GeV)");
    latex->DrawLatex(0.12,0.85,Form("%s (%s)",sSpecies.Data(), sCentLabels[iCent].Data()));


    for(Int_t in(0); in < listSyst->GetEntries(); ++in)
    {
      Int_t nnCol = gStyle->GetNumberOfColors();
      Int_t nPnt  = listSyst->GetEntries();
      Int_t idx = in * Float_t(nnCol-1) / (nPnt-1);

      TH1D* histSyst = (TH1D*) listSyst->At(in);
      if(!histSyst) { printf("no histSyst\n"); return;}
      StyleHist(histSyst,gStyle->GetColorPalette(idx),kOpenCircle);
      can->cd(1);
      histSyst->Draw("hist p same");
    }
    StyleHist(histSTD,kRed,kOpenSquare);
    // histSTD->SetFillColor(kRed);
    histSTD->Draw("same e1");

    can->cd(2);
    TH1* frame_2 = gPad->DrawFrame(histSTD->GetXaxis()->GetXmin(), 0.0,histSTD->GetXaxis()->GetXmax(), 0.2);
    frame_2->SetTitle("relative difference (RMS); p_{T} (GeV)");
    latex->DrawLatex(0.12,0.85,Form("%s (%s)",sSpecies.Data(), sCentLabels[iCent].Data()));
    leg->Draw("same");

    TH1D* histSigma = (TH1D*) histSTD->Clone(Form("sigma_cent%d",iCent));
    for(Int_t bin(1); bin < histSTD->GetNbinsX()+1; ++bin)
    {
      histSigma->SetBinContent(bin, histSTD->GetBinError(bin));
      histSigma->SetBinError(bin, 0.0);
    }

    // histSigma->SetMinimum(0.0);
    histSigma->Draw("same hist");

    can->SaveAs(Form("%s/%s_cent%d.pdf",sOutput.Data(),sSpecies.Data(), iCent), "pdf");

    fileOutput->cd();
    histSTD->Write(Form("mean_cent%d",iCent));
    histSigma->Write(Form("sigma_cent%d",iCent));
    listSyst->Write(Form("systs_cent%d",iCent),TObject::kSingleKey);

  }

}



// ==================================================================================================================
TH1D* GetRMS(TList* list)
{
  if(!list) {printf("no list\n"); return 0x0;}

  Int_t iEnt = list->GetEntries();

  TH1D* histRMS = (TH1D*) list->At(0)->Clone(Form("%s_std", list->At(0)->GetName()));
  if(!histRMS) {printf("no temp\n"); return 0x0;}

  histRMS->Reset();


  for(Int_t bin(1); bin < histRMS->GetNbinsX()+1; ++bin)
  {
    // mean
    Double_t dSum = 0.0;
    for(Int_t in(0); in < iEnt; ++in)
    {
      TH1D* temp = (TH1D*) list->At(in);
      Double_t dCon = temp->GetBinContent(bin);
      dSum += dCon;

      // printf(" %f ",dCon);
    }

    Double_t dMean = dSum / iEnt;
    // printf(" ||| %f ",dMean);

    // std
    Double_t dSumSq = 0.0;
    for(Int_t in(0); in < iEnt; ++in)
    {
      TH1D* temp = (TH1D*) list->At(in);
      Double_t dCon = temp->GetBinContent(bin);
      dSumSq += TMath::Power(dCon - dMean,2.0);
    }

    Double_t dSTD = dSumSq / (iEnt-1);
    // printf(" ||| %f \n",TMath::Sqrt(dSTD));

    // temp->SetBinConter(bin,TMath::Sqrt(dSTD));
    histRMS->SetBinContent(bin,dMean);
    histRMS->SetBinError(bin,TMath::Sqrt(dSTD));
  }

  return histRMS;
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

    if(bCor) { dErrRatio -= (2*dContNom*dErrDenom*dErrNom/TMath::Power(dContDenom,3)); }
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    ratio->SetBinContent(iBinNom,dContRatio);
    ratio->SetBinError(iBinNom,TMath::Sqrt(dErrRatio));
  }

  return ratio;
}
// ==================================================================================================================
TH1D* Diff(TH1* ratio)
{
  if(!ratio) { printf("Diff-ERROR: ratio histo does not exists\n"); return 0x0; }

  TH1D* hDiff = (TH1D*) ratio->Clone(Form("Diff_%s",ratio->GetName()));
  hDiff->Reset();

  for(Int_t bin(1); bin < hDiff->GetNbinsX()+1; ++bin)
  {
    Double_t dCont = ratio->GetBinContent(bin);
    Double_t dErr = ratio->GetBinError(bin);

    hDiff->SetBinContent(bin, TMath::Abs(dCont-1.0));
    hDiff->SetBinError(bin,dErr);
  }

  hDiff->SetMinimum(0.0);

  return hDiff;
}
// ==================================================================================================================
TH1D* BarlowTest(TH1* nom, TH1* denom, Bool_t bCor)
{
  if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

  Int_t binsNom = nom->GetNbinsX();
  Int_t binsDenom = denom->GetNbinsX();

  // if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  TH1D* hBarlow = (TH1D*) nom->Clone(Form("Barlow_%s_%s",nom->GetName(),denom->GetName()));
  hBarlow->Reset();

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


    Double_t dBarlowNom = TMath::Abs(dContDenom-dContNom);
    Double_t dBarlowDenom = 0.0;

    if(bCor) { dBarlowDenom = TMath::Abs(dErrDenom*dErrDenom - dErrNom*dErrNom); }
    else { dBarlowDenom = TMath::Abs(dErrDenom*dErrDenom + dErrNom*dErrNom); }
    // Double_t dBarlowDenom = dErrDenom*dErrDenom - dErrNom*dErrNom;

    if(dBarlowDenom <= 0.0) continue;

    hBarlow->SetBinContent(iBinDenom,dBarlowNom/TMath::Sqrt(dBarlowDenom));
  }

  hBarlow->SetMinimum(0.0);

  return hBarlow;
}
// ==================================================================================================================
TH1D* ApplyBarlow(TH1* diff, TH1* barlow, Double_t dCut)
{
  if(!diff) { printf("ERROR-ApplyBarlow : no diff\n"); return 0x0; }
  if(!barlow) { printf("ERROR-ApplyBarlow : no barlow\n"); return 0x0; }

  TH1D* after = (TH1D*) diff->Clone(Form("%s_Barlowed",diff->GetName()));
  if(dCut < 0.0) { return after; }

  for(Int_t iBin(1); iBin < after->GetNbinsX()+1; ++iBin)
  {
    if(barlow->GetBinContent(iBin) <= dCut)
    {
      after->SetBinError(iBin, 999.9);
    }
  }

  return after;
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
// ==================================================================================================================
void SetCustomPalette()
{
  // Setting custom (ala ROOT6) color palette
  // See https://root.cern.ch/doc/master/TColor_8cxx_source.html#l02400
  // for definition of color setting (array bellow) for ROOT 6 defined palettes

  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000,	0.6250, 0.7500, 0.8750, 1.0000};

  // Rain Bow
  Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
  Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
  Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};


  // Bird
  //case 57:
  // Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  // Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  // Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};

  // Blue Green Yellow
  // //case 71:
  // Double_t red[9]   = { 22./255., 19./255.,  19./255.,  25./255.,  35./255.,  53./255.,  88./255., 139./255., 210./255.};
  // Double_t green[9] = {  0./255., 32./255.,  69./255., 108./255., 135./255., 159./255., 183./255., 198./255., 215./255.};
  // Double_t blue[9]  = { 77./255., 96./255., 110./255., 116./255., 110./255., 100./255.,  90./255.,  78./255.,  70./255.};

  // Viridis
  // case 112:
  // Double_t red[9]   = { 26./255., 51./255.,  43./255.,  33./255.,  28./255.,  35./255.,  74./255., 144./255., 246./255.};
  // Double_t green[9] = {  9./255., 24./255.,  55./255.,  87./255., 118./255., 150./255., 180./255., 200./255., 222./255.};
  // Double_t blue[9]  = { 30./255., 96./255., 112./255., 114./255., 112./255., 101./255.,  72./255.,  35./255.,   0./255.};

  // Cividis
  // case 113:
  // Double_t red[9]   = {  0./255.,   5./255.,  65./255.,  97./255., 124./255., 156./255., 189./255., 224./255., 255./255.};
  // Double_t green[9] = { 32./255.,  54./255.,  77./255., 100./255., 123./255., 148./255., 175./255., 203./255., 234./255.};
  // Double_t blue[9]  = { 77./255., 110./255., 107./255., 111./255., 120./255., 119./255., 111./255.,  94./255.,  70./255.};

  Int_t pal = TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, 1);
  const Int_t nCol = 255;
  Int_t colors[nCol];
  for (int i=0; i<nCol; i++) colors[i] = pal+i;

  gStyle->SetPalette(nCol,colors);
}
