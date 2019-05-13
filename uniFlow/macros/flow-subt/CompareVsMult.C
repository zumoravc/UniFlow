#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"

Color_t gColors[] = { kRed, kGreen+2, kBlue, kMagenta+2, kBlack};
Color_t gMarkerStyles[] = { kFullCircle, kOpenCircle, kOpenSquare,kOpenSquare,kOpenSquare };

const char* gsInputFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/high_pt_range/output_weighted/";
// const char* gsInputFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_noweight/";
// const char* gsInputFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_sub/";
// const char* gsInputFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_sub_norm/";
// const char* gsInputFolder = "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_sub_norm_weighted/";

// const char* gsFileName = "Subtracted_Charged.root";
const char* gsFileName = "Raw_Charged.root";

const char* gsOutputFolder = gsInputFolder;

// const char* gsFileList[] = {
//   "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_weighted/Subtracted_Charged.root",
//   "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_noweight/Subtracted_Charged.root",
//   "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_sub/Subtracted_Charged.root",
//   "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_sub_norm/Subtracted_Charged.root",
//   "/Users/vpacik/NBI/Flow/uniFlow/results/flowsub/pPb-run2-gap0/output_sub_norm_weighted/Subtracted_Charged.root"
// };

const char* gsLegendLabels[] = {"0-20%","20-40%","40-60%","60-100%"};

const char* gsHistoName = "hFlow2_Charged_harm2_gap00";
// Double_t gdMinimum = -0.15;
// Double_t gdMaximum = 0.4;

Double_t gdMinimum = 0.;
Double_t gdMaximum = 1.;


void CompareVsMult()
{

  TString sFileName = Form("%s/%s",gsInputFolder,gsFileName);
  Int_t iNumCent = 4;

  // output files
  TFile* fInputFile = TFile::Open(sFileName.Data(),"READ");
  if(!fInputFile) { printf("ERROR: Input file not found.\n"); return; }

  fInputFile->ls();

  TLegend* leg = new TLegend(0.15,0.6,0.4,0.88);
  leg->SetBorderSize(0);


  TCanvas* can = new TCanvas();

  for(Int_t iCent(0); iCent < iNumCent; ++iCent)
  {
    TH1D* hTemp = (TH1D*) fInputFile->Get(Form("%s_cent%d",gsHistoName,iCent));
    if(!hTemp) { printf("ERROR: Temp histo not found!\n"); continue; }

    leg->AddEntry(hTemp,gsLegendLabels[iCent],"pl");
    hTemp->SetStats(0);
    hTemp->SetTitle("Centrality dependence");
    hTemp->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    hTemp->GetYaxis()->SetTitle("v_{2} (unsub)");
    hTemp->SetMaximum(gdMaximum);
    hTemp->SetMinimum(gdMinimum);
    hTemp->SetLineColor(gColors[iCent]);
    hTemp->SetMarkerColor(gColors[iCent]);
    hTemp->SetMarkerStyle(kFullSquare);
    can->cd();
    hTemp->Draw("same");


  }
  leg->Draw();
  can->SaveAs(Form("%s/cent_dependence.pdf",gsOutputFolder),"pdf");


  return;
}
