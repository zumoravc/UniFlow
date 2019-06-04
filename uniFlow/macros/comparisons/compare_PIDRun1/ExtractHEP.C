


#include "TFile.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"


TH1F* TranslateHist(TH1F* hist = 0x0, TH1F* histErr = 0x0);

void ExtractHEP()
{
  const char* sInputHEP = "/Users/vpacik/NBI/Flow/data/HEPdata/PLB726(2013)/HEPData-ins1242302-v1-root.root";
  const char* sOutput = "/Users/vpacik/NBI/Flow/data/HEPdata/PLB726(2013)/extracted-v2SP.root";

  TFile* fileHEP = new TFile(sInputHEP,"READ");
  if(!fileHEP) { printf("HEP file not found.\n"); return; }

  TFile* fileOutput = new TFile(sOutput,"RECREATE");
  if(!fileOutput) { printf("Output file not found.\n"); return; }

  fileHEP->ls();

  const Int_t iNumCent = 4;
  Int_t iTable[iNumCent] = {5,9,13,17};

  const Int_t iNumSpecies = 4;
  TString sSpecies[iNumSpecies] = {"Charged","Pion","Kaon","Proton"};

  TH1F* hist = 0x0;
  TH1F* histErr = 0x0;
  TH1F* histDone = 0x0;
  TGraphAsymmErrors* graph = 0x0;

  TCanvas* cfig = new TCanvas("cfig", "Alice Figure Template", 800, 600);
  TH1* h = cfig->DrawFrame(0,0,6,0.4);

  for(Int_t iCent(0); iCent < iNumCent; iCent++)
  {
    for(Int_t iSpecies(0); iSpecies < iNumSpecies; iSpecies++)
    {
      printf("%d\n",iTable[iCent]+iSpecies);

      fileHEP->cd(Form("Table %d",iTable[iCent]+iSpecies));
      gDirectory->ls();

      graph = (TGraphAsymmErrors*) gDirectory->Get("Graph1D_y2");
      if(!graph) { printf("No graph found\n"); return; }
      hist = (TH1F*) gDirectory->Get("Hist1D_y2");
      if(!hist) { printf("No hist found\n"); return; }
      histErr = (TH1F*) gDirectory->Get("Hist1D_y2_e1");
      if(!histErr) { printf("No histErr found\n"); return; }

      // hist->Draw("same hist");
      // histErr->Draw("same p");
      graph->Draw("same p");

      histDone = TranslateHist(hist,histErr);
      if(!histDone) { printf("hist done not done properly\n"); return; }
      histDone->SetName(Form("hFlow2_%s_harm2_gap08_cent%d",sSpecies[iSpecies].Data(),iCent));

      histDone->Draw("same p2");

      fileOutput->cd();
      histDone->Write();
      // histDone->Write(Form("hFlow2_%s_harm2_gap08_cent%d",sSpecies[iSpecies].Data(),iCent));
    }
  }

  histDone->ls();



  // graph->Draw("same p2");



  return;
}

TH1F* TranslateHist(TH1F* hist, TH1F* histErr)
{
  if(!hist || !histErr) { printf("At least one of input not found\n"); return 0x0; }

  Int_t iNumBins = hist->GetNbinsX();
  Int_t iNumBinsErr = histErr->GetNbinsX();
  // printf("%d | %d\n",hist->GetNbinsX(),histErr->GetNbinsX());
  printf("%d | %d\n",iNumBins,iNumBinsErr);

  if(iNumBins != iNumBinsErr) { printf("Not the same number of bins\n"); return 0x0; }

  TH1F* histDone = (TH1F*) hist->Clone("histDone");
  for(Int_t iBin(1); iBin < iNumBins+1; iBin++)
  {
    histDone->SetBinError(iBin,histErr->GetBinError(iBin));
    // hist->SetBinError(iBin,histErr->GetBinError(iBin));
  }

  return histDone;
  // return hist;
}
