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

// colors
Color_t colors[] = {kRed, kGreen+2, kBlue, kMagenta+1};
Style_t markers[] = {kOpenSquare, kOpenCircle, kOpenDiamond};

void Plot_vectorMethods()
{
  const Int_t iNumMethods = 3;
  TString sMethods[iNumMethods] = {"GF_eventweighted", "GF_noneventweighted", "SP_nonscaled_noneventweighted"};
  TString sMethodLabel[iNumMethods] = {"Generic framework (GF)", "GF w/o event weights", "GF w/o event weigts & not normalized vectors"};
  TString sSpecies = "Charged";
  TString sListName = "list_SubtPPb_vn";
  TString sHistName = "hSubPPb_vn_cent";
  TString sLegendHeader = "pPb - pPb (60-100%)";

  TString sOutFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gap08/output_vn/vectorMethod/";

  Int_t iNumCent = 4;
  TString sCentLabel[] = {"0-20%", "20-40%", "40-60%", "60-100%"};

  gSystem->mkdir(sOutFolder.Data(), kTRUE);
  for(Int_t cent(0); cent < iNumCent; ++cent)
  {

    TLegend* leg = new TLegend(0.12,0.65,0.65,0.88);
    leg->SetHeader(sLegendHeader.Data());
    leg->SetBorderSize(0.);
    leg->SetFillColor(0);

    TCanvas* can = new TCanvas("can","can",400,400);
    can->cd();
    TH1* frame_can = (TH1*) gPad->DrawFrame(0.,0.,10.,0.3);
    frame_can->SetTitle(Form("v_{2}{2}^{sub} (%s V0A); p_{T} (GeV/c)",sCentLabel[cent].Data()));

    for(Int_t method(0); method < iNumMethods; ++method)
    {
      TString sFile = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run3-gap08/output_vn/" + sMethods[method] + "/pPb08_pp08_subt/" + sSpecies + "/Subt_results.root";

      TFile* fileIn = OpenFile(sFile.Data()); if(!fileIn) { return; }
      fileIn->ls();

      TList* listIn = (TList*) fileIn->Get(sListName.Data()); if(!listIn) { return; }
      listIn->ls();

      TH1D* hist = (TH1D*) listIn->FindObject(Form("%s%d",sHistName.Data(),cent)); if(!hist) { return; }
      StyleHist(hist,colors[method],markers[method]);

      leg->AddEntry(hist,sMethodLabel[method].Data(),"p");

      can->cd();
      leg->Draw();
      hist->Draw(" e1 same");

    }

    // can->SaveAs(Form("%s%s_cent%d.pdf",sOutputFolder.Data(),sListName.Data(),cent),"pdf");
    can->SaveAs(Form("%s%s_cent%d.pdf",sOutFolder.Data(),sListName.Data(),cent),"pdf");
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
