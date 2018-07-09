#include "TFile.h"
#include "TH1.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"

void StyleHist(TH1* hist, Color_t color = kRed, Style_t markerStyle = kOpenCircle);


Color_t colorsMethods[] = {kRed, kGreen+2, kBlue+2};
Color_t markersMethods[] = {kOpenSquare, kOpenCircle, kOpenDiamond};
TString sMethods[] = {"GF_eventweighted","GF_noneventweighted","SP_nonscaled_noneventweighted"};
TString sMethodLabels[] = {"scaled & event weights (A)","scaled & W/O event weights (B)"," NOT scaled & W/O event weigths (C)"};


Color_t colorsBases[] = {kGreen+2, kBlue+2, kRed};
Color_t markersBases[] = {kFullCircle, kFullCircle, kFullCircle};
TString sBaseLines[] = {"pPb-pp-cent", "pPb-pp-int","pPb-pPb-perp"};
TString sBaseLineLabels[] = {Form("pp (%s)","same cent."), "pp (MB, 0-100%)","pPb (60-100%)"};

TString sGaps[] = {"gap04"};
Double_t dGapValues[] = {0.4};
TString sSpecies_list[] = {"Charged"};

const Int_t iNumCent = 4;
TString sCentLabels[iNumCent] = {"0-20%", "20-40%", "40-60%", "60-100%"};

TString sPathTop = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/etegap-dependence/output-comp/";

Double_t dPtmin = 0.0;
Double_t dPtmax = 10.0;
Double_t dYmin = -0.005;
Double_t dYmax = 0.7;

void PlotBaselines()
{
  Int_t iNumMethods = sizeof(sMethods)/sizeof(sMethods[0]);
  Int_t iNumBases = sizeof(sBaseLines)/sizeof(sBaseLines[0]);

  Int_t iGap = 0;
  Int_t iSpecies = 0;


  TString sSpecies = sSpecies_list[iSpecies];
  TString sGap = sGaps[iGap];
  Double_t dGapValue = dGapValues[iGap];

  for(Int_t iMethod(0); iMethod < iNumMethods; ++iMethod)
  {
    TString sMethod = sMethods[iMethod];
    TString sMethodLabel = sMethodLabels[iMethod];

    for(Int_t iCent = 0; iCent < iNumCent; ++iCent)
    {
      TString sCentLabel = sCentLabels[iCent];

      TCanvas* canBases = new TCanvas("canBases","canBases", 400,400);
      canBases->cd(1);
      TH1* frame_canBases = (TH1*) gPad->DrawFrame(dPtmin,dYmin,dPtmax,dYmax);
      frame_canBases->SetTitle(Form("%s (%s); p_{T} (GeV/c); v^{sub}_{2}{2,|#Delta#eta| < %.1f}",sSpecies.Data(), sCentLabel.Data(),dGapValue));

      TLegend* legMethods = new TLegend(0.12,0.55,0.4,0.85);
      legMethods->SetFillColorAlpha(0,0);
      legMethods->SetBorderSize(0);
      legMethods->SetHeader(Form("%s",sMethodLabel.Data()));

      for(Int_t iBase = 0; iBase < iNumBases; ++iBase)
      {
        TString sBase = sBaseLines[iBase];

        TString sInFile = sPathTop + sGap + "/" + sMethod + "/" + sBase + "/Subt_results.root";

        TFile* fileIn = TFile::Open(sInFile.Data(),"READ"); if(!fileIn) { printf("No fileIn\n"); return; }

        if(iBase == 0)
        {
          // un-subtracted
          TH1D* raw = (TH1D*) fileIn->Get(Form("hVn_Raw_cent%d",iCent));
          if(!raw) { printf("No hist\n"); fileIn->ls(); return; }
          StyleHist(raw, kBlack, kOpenSquare);
          legMethods->AddEntry(raw,Form("Unsubt p-Pb"),"p");

          canBases->cd(1);
          raw->DrawCopy("same");
        }

        TH1D* hist = (TH1D*) fileIn->Get(Form("hVn_Subt_cent%d",iCent));
        if(!hist) { printf("No hist\n"); fileIn->ls(); return; }
        StyleHist(hist, colorsBases[iBase],markersBases[iBase]);
        legMethods->AddEntry(hist,sBaseLineLabels[iBase].Data(),"p");

        canBases->cd(1);
        hist->DrawCopy("same");

      }
      canBases->cd(1);
      legMethods->Draw();

      TString sOutput = Form("%s/comp-bases/%s/",sPathTop.Data(),sGap.Data());
      gSystem->mkdir(sOutput.Data(),kTRUE);
      canBases->SaveAs(Form("%s/vn_subt_%s_cent%d.pdf",sOutput.Data(),sMethod.Data(),iCent),"pdf");
    }
  }

  return;
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
