#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TSystem.h"
#include "TMath.h"



void PlotVnPtSubt()
{
  TString sTag = "weighted";
  // TString sTag = "noweight";
  // TString sTag = "sub";
  // TString sTag = "sub_norm";
  // TString sTag = "sub_norm_weighted";


  TString sGap = "00";
  TString sInputFile = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap00/output_%s/flow-subtraction/Subtracted.root",sTag.Data());
  TString sOutputDir = Form("/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap00/output_%s/flow-subtraction/",sTag.Data());
  Int_t iNumCent = 4;

  Int_t iNumSpecies = 7;
  TString sSpecies[] = {"Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"};
  Color_t gColors[] = { kBlack, kRed, kGreen+2, kCyan+1 , kBlue, kOrange+1, kMagenta};
  Color_t gMarkerStyles[] = { kFullCircle, kFullStar, kFullTriangleUp, kFullSquare,kFullTriangleDown, kFullDiamond, kFullCross };

  Double_t dXmin = 0.2;
  Double_t dXmax = 7.;
  Double_t dYmin = -0.1;
  Double_t dYmax = 0.3;

  TFile* fInputFile = TFile::Open(sInputFile.Data(),"READ");
  if(!fInputFile) { printf("ERROR: Input file '%s' not found!\n",sInputFile.Data()); return; }

  for(Int_t iCent(0); iCent < iNumCent; ++iCent)
  {
    TLegend* leg = new TLegend(0.12,0.6,0.4,0.89);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0.,0.);
    TCanvas* canvas = new TCanvas("canvas","canvas");
    canvas->cd();
    TH1F* hFrame = canvas->DrawFrame(dXmin,dYmin,dXmax,dYmax,"; p_{T} (GeV/c); v_{2}^{sub} {|#Delta#eta| > 0}");

    for(Int_t iSpecies(0); iSpecies < iNumSpecies; ++iSpecies)
    {
      TString sHistoName = Form("hFlow2_%s_harm2_gap%s",sSpecies[iSpecies].Data(),sGap.Data());
      if( sSpecies[iSpecies].EqualTo("Charged") ||
        sSpecies[iSpecies].EqualTo("Pion") ||
        sSpecies[iSpecies].EqualTo("Kaon") ||
        sSpecies[iSpecies].EqualTo("Proton") )
      {
        sHistoName.Append(Form("_cent%d_subt",iCent));
      }
      else
      {
        sHistoName.Append(Form("_mult%d_subt",iCent));
      }

      TH1D* hTemp = (TH1D*) fInputFile->Get(sHistoName.Data());
      if(!hTemp) { printf("WARNING: Histo '%s' not found\n",sHistoName.Data()); continue; }
      hTemp->SetLineColor(gColors[iSpecies]);
      hTemp->SetMarkerColor(gColors[iSpecies]);
      hTemp->SetMarkerStyle(gMarkerStyles[iSpecies]);
      hTemp->Draw("same");
      leg->AddEntry(hTemp,sSpecies[iSpecies].Data(),"pl");
    }

    leg->Draw();
    canvas->SaveAs(Form("%s/VnPtSub_cent%d.pdf",sOutputDir.Data(),iCent),"pdf");
  }






  return;
}
