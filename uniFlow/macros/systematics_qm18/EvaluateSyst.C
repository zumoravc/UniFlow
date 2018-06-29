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
#include "TColor.h"
#include "TSystem.h"
#include "TStyle.h"

// void Eval_single(TString sSpecies, Double_t* dDiffs, Bool_t bPtDep);
void SetCustomPalette();


TString sCentLabel[] = {"0-10%","10-20%","20-40%","40-60%"};
Int_t iNumCent = 4;

TString sInputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/out-binning-3/";
TString sOutputPath = "/Users/vpacik/NBI/Flow/uniFlow/results/qm-run/syst/out-binning-3/final/";

TH1D* histTotal = 0x0;
Int_t iNumBins = 0;


Color_t col[] = {kRed, kGreen+2, kBlue, kMagenta-1, kOrange+2, kSpring+1};

void EvaluateSyst()
{

  // TString sSpecies = "Charged"; TString sSpeciesFB = "Charged";  Bool_t bFB = kTRUE; Bool_t bFitting = kFALSE; Bool_t bUsePalette = kTRUE;
  // TString sLabels[] = {"PV","TPCcls","FB"};  Double_t dDiffs[4][2] = {{2,1},{2,1},{2,1},{2,1}};

  // TString sSpecies = "Pion"; TString sSpeciesFB = "Pion";  Bool_t bFB = kTRUE; Bool_t bFitting = kFALSE;Bool_t bUsePalette = kTRUE;
  // TString sLabels[] = {"PV","TPCcls","PID (Bayes)","PID(nsgima)"};  Double_t dDiffs[4][4] = {{2,1,0,2},{2,1,0,2},{2,1,0,2},{2,1,0,2}};

  // TString sSpecies = "Kaon"; TString sSpeciesFB = "Pion";  Bool_t bFB = kTRUE; Bool_t bFitting = kFALSE;Bool_t bUsePalette = kTRUE;
  // TString sLabels[] = {"PV","TPCcls","PID (Bayes)","PID(nsgima)"};  Double_t dDiffs[4][4] = {{2,2,2,2},{2,2,2,2},{4,2,0,4},{4,2,0,5}};

  // TString sSpecies = "Proton"; TString sSpeciesFB = "Pion";  Bool_t bFB = kTRUE; Bool_t bFitting = kFALSE;Bool_t bUsePalette = kTRUE;
  // TString sLabels[] = {"PV","TPCcls","PID (Bayes)","PID(nsgima)"};  Double_t dDiffs[4][4] = {{2,2,0,4},{2,2,0,4},{4,2,0,4},{4,2,0,6}};

  // TString sSpecies = "K0s"; TString sSpeciesFB = ""; Bool_t bFB = kFALSE; Bool_t bFitting = kTRUE; Bool_t bUsePalette = kTRUE;
  // TString sLabels[] = {"PV","TPCcls","FB","DCA","Decay radius","CPA","PID (nsigma)",};  Double_t dDiffs[4][7] = {{4,2,7,1,2,2,2},{4,2,7,1,2,2,2},{4,2,7,2,3,3,2},{6,2,7,4,5,6,5}};

  // TString sSpecies = "Lambda"; TString sSpeciesFB = ""; Bool_t bFB = kFALSE; Bool_t bFitting = kTRUE; Bool_t bUsePalette = kTRUE;
  // TString sLabels[] = {"PV","TPCcls","FB","DCA","Decay radius","CPA","PID (nsigma)",};  Double_t dDiffs[4][7] = {{4,2,7,1,2,2,2},{4,2,7,2,2,3,2},{4,2,7,2,3,3,2},{6,2,7,4,4,6,5}};

  TString sSpecies = "Phi"; TString sSpeciesFB = ""; Bool_t bFB = kFALSE; Bool_t bFitting = kTRUE; Bool_t bUsePalette = kTRUE;
  TString sLabels[] = {"PV","TPCcls","FB","PID (Bayes)","PID (nsigma)",};  Double_t dDiffs[4][5] = {{4,2,7,2,0},{4,2,7,2,0},{4,2,7,4,0},{6,2,7,5,0}};

  Int_t iNumDiffs = sizeof(dDiffs[0]) / sizeof(dDiffs[0][0]);
  // printf("%d",iNumDiffs);

  // ======================

  gSystem->mkdir(sOutputPath.Data(),1);
  TFile* fileOutput = TFile::Open(Form("%s/syst_%s.root",sOutputPath.Data(),sSpecies.Data()),"RECREATE");
  if(!fileOutput) return;

  for(Int_t iCent = 0; iCent < iNumCent; ++iCent)
  {

    TLegend* leg = new TLegend(0.5,0.45,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0,0);
    leg->SetHeader(Form("%s (%s)",sSpecies.Data(), sCentLabel[iCent].Data()));

    TList* listSyst = new TList();

    for(Int_t iDiff = 0; iDiff < iNumDiffs; ++iDiff)
    {
      Double_t dValue = dDiffs[iCent][iDiff] / 100.0;

      // extract binning for given species (only once)
      if(iDiff == 0)
      {
        TString sInputFile = sInputPath + "tracking/FB.root";

        TFile* fileInput = TFile::Open(sInputFile.Data(),"READ");
        if(!fileInput) { return; }

        TList* list = (TList*) fileInput->Get(Form("%s_cent%d",sSpecies.Data(),iCent));
        if(!list) { printf("no list\n"); return; }
        // list->ls();
        TH1D* hist = (TH1D*) list->FindObject("histBase");
        if(!hist) { printf("no hist\n"); return; }

        histTotal = (TH1D*) hist->Clone(Form("total_%s_cent%d",sSpecies.Data(), iCent));
        if(!histTotal) { printf("no hist syst\n"); return; }

        iNumBins = histTotal->GetNbinsX();
        histTotal->Reset();
      }

      TH1D* histSyst = (TH1D*) histTotal->Clone(Form("syst_%s_%s_cent%d",sSpecies.Data(), sLabels[iDiff].Data(), iCent));
      histSyst->Reset();
      leg->AddEntry(histSyst,sLabels[iDiff].Data(),"l");

      for(Int_t bin(1); bin < iNumBins+1; ++bin)
      {
        Double_t dOld = histTotal->GetBinContent(bin);
        Double_t dNew = dValue*dValue + dOld*dOld;

        histTotal->SetBinContent(bin, TMath::Sqrt(dNew));
        histSyst->SetBinContent(bin, dValue);

      }
      listSyst->Add(histSyst);
    }

    // pT dependent

    // FB
    if(bFB)
    {
      TFile* fileFB = TFile::Open(Form("%s/tracking/FB.root",sInputPath.Data()),"READ");
      if(!fileFB) { return; }
      TList* listFB = (TList*) fileFB->Get(Form("%s_cent%d",sSpeciesFB.Data(),iCent));
      TF1* fitRatioDiff = (TF1*)listFB->FindObject("fitRatioDiff");
      if(!fitRatioDiff) {printf("no fb fit"); return; }

      TH1D* histFB = (TH1D*) histTotal->Clone(Form("syst_%s_%s_cent%d",sSpecies.Data(), "FB", iCent));
      histFB->Reset();
      leg->AddEntry(histFB,"FB","l");

      for(Int_t bin(1); bin < iNumBins+1; ++bin)
      {
        Double_t dBinCenter = histFB->GetBinCenter(bin);

        Double_t dUnc = 0.0;
        if(dBinCenter <= fitRatioDiff->GetXmax())
        {
          dUnc = fitRatioDiff->Eval(dBinCenter);
        }

        histFB->SetBinContent(bin,dUnc);

        Double_t dOld = histTotal->GetBinContent(bin);
        Double_t dNew = dUnc*dUnc + dOld*dOld;
        histTotal->SetBinContent(bin, TMath::Sqrt(dNew));
      }

      listSyst->Add(histFB);
    }

    if(bFitting)
    {
      TFile* fileFit = TFile::Open(Form("%s/fitting/fitting_%s.root",sInputPath.Data(),sSpecies.Data()),"READ");
      if(!fileFit) { printf("no fileFit\n"); return; }
      TH1D* histFit = (TH1D*) fileFit->Get(Form("sigma_cent%d",iCent));
      if(iCent == 3) { histFit = (TH1D*) fileFit->Get(Form("sigma_cent%d",2)); }
      if(!histFit) { printf("no histFit\n"); fileFit->ls(); return; }


      leg->AddEntry(histFit,"fitting","l");
      listSyst->Add(histFit);

      for(Int_t bin(1); bin < iNumBins+1; ++bin)
      {
        Double_t dUnc = histFit->GetBinContent(bin);

        Double_t dOld = histTotal->GetBinContent(bin);
        Double_t dNew = dUnc*dUnc + dOld*dOld;
        histTotal->SetBinContent(bin, TMath::Sqrt(dNew));
      }

    }

    // here all contrib in list;

    // saving
    fileOutput->cd();
    listSyst->Write(Form("cent%d",iCent),TObject::kSingleKey);
    histTotal->Write(Form("syst_cent%d",iCent));


    // now drawing

    SetCustomPalette();
    Int_t nnCol = gStyle->GetNumberOfColors();
    Int_t nPnt  = listSyst->GetEntries();


    TCanvas* can = new TCanvas("can","can",600,600);
    can->cd();
    // TH1* frame_can = (TH1*) gPad->DrawFrame(histTotal->GetXaxis()->GetXmin(), 0.0, histTotal->GetXaxis()->GetXmax(), histTotal->GetMaximum()+0.05);
    TH1* frame_can = (TH1*) gPad->DrawFrame(histTotal->GetXaxis()->GetXmin(), 0.0, histTotal->GetXaxis()->GetXmax(), 0.3);
    frame_can->SetTitle("; p_{T} (GeV/c); Relative uncertainy");

    for(Int_t i(0); i < nPnt; ++i)
    {
      TH1D* hist = (TH1D*) listSyst->At(i);
      if(!hist) { printf("no form list \n"); return; }

      Int_t idx = i * Float_t(nnCol-1) / (nPnt-1);

      can->cd();
      if(bUsePalette){ hist->SetLineColor(gStyle->GetColorPalette(idx)); }
      else {  hist->SetLineColor(col[i]); }


      hist->Draw("hist same");
    }

    can->cd();
    histTotal->SetLineColor(kRed);
    histTotal->SetLineWidth(2);



    histTotal->DrawCopy("hist same");
    leg->AddEntry(histTotal, "Total","l");
    leg->Draw();

    can->SaveAs(Form("%s/syst_%s_cent%d.pdf",sOutputPath.Data(), sSpecies.Data(), iCent),"pdf");
  }
}



// void Eval_single(TString sSpecies, Double_t* dDiffs, Bool_t bPtDep)
// {
//
//   TFile* fileInput = TFile::Open(sInputFile.Data(),"READ");
//   if(!fileInput) { return; }
//
//   gSystem->mkdir(sOutputPath.Data(),1);
//   TFile* fileOutput = TFile::Open(Form("%s/syst_%s.root",sOutputPath.Data(),sSpecies.Data()),"RECREATE");
//   if(!fileOutput) { return; }
//
//   // fileInput->ls();
//
//   for(Int_t iCent(0); iCent < iNumCent; ++iCent)
//   {
//     TList* list = (TList*) fileInput->Get(Form("%s_cent%d",sSpecies.Data(),iCent));
//     if(!list) { printf("no list\n"); return; }
//     TH1D* hist = (TH1D*) list->FindObject("histBase");
//     if(!hist) { printf("no hist\n"); return; }
//
//
//     TH1D* histSyst = (TH1D*) hist->Clone(Form("syst_%s_cent%d",sSpecies.Data(), iCent));
//     if(!histSyst) { printf("no hist syst\n"); return; }
//
//     TF1* fitDiff = 0x0;
//     if(bPtDep) { fitDiff = (TF1*) list->FindObject("fitRatioDiff"); if(!fitDiff) { printf("no fit diff\n"); return; } }
//
//     histSyst->Reset();
//
//     for(Int_t iBin(1); iBin<histSyst->GetNbinsX()+1; ++iBin)
//     {
//       Double_t dUnc = dDiffs[iCent]*dDiffs[iCent];
//
//       if(bPtDep)
//       {
//         Double_t dBinCenter = histSyst->GetBinCenter(iBin);
//
//         if(dBinCenter <= fitDiff->GetXmax())
//         {
//           Double_t dDiff = fitDiff->Eval(dBinCenter);
//           dUnc += dDiff*dDiff;
//         }
//       }
//
//       histSyst->SetBinContent(iBin, TMath::Sqrt(dUnc));
//     }
//
//     histSyst->Draw();
//     fileOutput->cd();
//     histSyst->Write();
//
//   }
//
//
//   return;
//
// }
void SetCustomPalette()
{
  // Setting custom (ala ROOT6) color palette
  // See https://root.cern.ch/doc/master/TColor_8cxx_source.html#l02400
  // for definition of color setting (array bellow) for ROOT 6 defined palettes

  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000,	0.6250, 0.7500, 0.8750, 1.0000};

  // // Rain Bow
  Double_t red[9]   = {  0./255.,   5./255.,  15./255.,  35./255., 102./255., 196./255., 208./255., 199./255., 110./255.};
  Double_t green[9] = {  0./255.,  48./255., 124./255., 192./255., 206./255., 226./255.,  97./255.,  16./255.,   0./255.};
  Double_t blue[9]  = { 99./255., 142./255., 198./255., 201./255.,  90./255.,  22./255.,  13./255.,   8./255.,   2./255.};

  // Visible Spectrum
  // Double_t red[9]   = { 18./255.,  72./255.,   5./255.,  23./255.,  29./255., 201./255., 200./255., 98./255., 29./255.};
  // Double_t green[9] = {  0./255.,   0./255.,  43./255., 167./255., 211./255., 117./255.,   0./255.,  0./255.,  0./255.};
  // Double_t blue[9]  = { 51./255., 203./255., 177./255.,  26./255.,  10./255.,   9./255.,   8./255.,  3./255.,  0./255.};


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
