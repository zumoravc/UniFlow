/*
* Macro for plotting the comparison between FAST, CENTwSDD and CENTwoSDD trigger clusters
* QA and flow results.
* It contains comparison between flow from POIs in positive and negative eta as requested by Flow PAG
* Vojtech Pacik (vojtech.pacik@cern.ch) - 30 May 2017
*/

#include <vector>
#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"


const Short_t fiNumFiles = 3;
TString fsOutputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_Recon_full_16q/comparison";
TString fsOutputFileFormat = "pdf";
// convention: everything in arrays: / 0: FAST / 1: CENTwSDD / 2: CENTwoSDD
TString fsFileName[fiNumFiles] = {"/Users/vpacik/NBI/Flow/results/uniFlow_Recon_full_16q/FAST/AnalysisResults.root", "/Users/vpacik/NBI/Flow/results/uniFlow_Recon_full_16q/CENT_wSDD/AnalysisResults.root", "/Users/vpacik/NBI/Flow/results/uniFlow_Recon_full_16q/CENT_woSDD/AnalysisResults.root"};
TString fsFileNameFlow[fiNumFiles] = {"/Users/vpacik/NBI/Flow/results/uniFlow_Recon_full_16q/FAST/results/UniFlow.root", "/Users/vpacik/NBI/Flow/results/uniFlow_Recon_full_16q/CENT_wSDD/results/UniFlow.root", "/Users/vpacik/NBI/Flow/results/uniFlow_Recon_full_16q/CENT_woSDD/results/UniFlow.root"};
TString fsTriggerName[fiNumFiles] = {"FAST","CENT_wSDD","CENT_woSDD"};

// defining colors
Color_t fColors[fiNumFiles] = {kRed, kBlue, kGreen-2};
Short_t fMarkers[2*fiNumFiles] = {kFullSquare, kFullCircle,  kFullTriangleUp, kOpenSquare, kOpenCircle,  kOpenTriangleUp};

TCanvas* PlotList(const TList* list = 0x0, const Short_t iPlotRatio = 0, const char* sDrawOpt = "hist");
TCanvas* PlotListPosNeg(const TList* list = 0x0, const Short_t iPlotRatio = 0, const char* sDrawOpt = "hist");
void PlotComparison(TString histName = "", TString histListName = "", const char* histClass = "TH1D");


void PlotComparisonCentFast()
{
  // Lists to be processed
  TList* list_fhEventCentrality = new TList();
  TList* list_fhQAV0sInvMassK0s = new TList();
  TList* list_fhQAV0sInvMassLambda = new TList();
  TList* list_fhPhiInvMass = new TList();
  TList* list_fhRefsMult = new TList();
  TList* list_fhRefsPt = new TList();
  TList* list_fhRefsPhi = new TList();
  TList* list_fhRefsEta = new TList();
  TList* list_hFlow2_Refs_harm2_gap08 = new TList();
  TList* list_hFlow2_Charged_harm2_gap08_cent0 = new TList();
  TList* list_hFlow2_Charged_harm2_gap08_cent1 = new TList();
  TList* list_hFlow2_Charged_harm2_gap08_cent2 = new TList();
  TList* list_hFlow2_Charged_harm2_gap08_cent3 = new TList();
  TList* list_hFlow2_Charged_harm2_gap08_cent4 = new TList();


  // placeholders
  TFile* file = 0x0;
  TList* list = 0x0;
  TH1* hist = 0x0;

  // openning files & loading various histos
  for(Short_t iFile(0); iFile < fiNumFiles; iFile++)
  {
    // QA histograms
    file = new TFile(fsFileName[iFile].Data(),"READ");
    if(!file) { printf("Input file '%s' not loaded\n",fsFileName[iFile].Data()); return; }
    file->cd("UniFlow");

    // loading histos

    // Events
    list = (TList*) gDirectory->Get("QA_Events_UniFlow");
    if(!list) { printf("List not loaded\n"); return; }

    hist = (TH1D*) list->FindObject("fhEventCentrality");
    hist->Scale(1/hist->GetEntries());
    if(hist) list_fhEventCentrality->Add(hist);

    // // Charged
    list = (TList*) gDirectory->Get("QA_Charged_UniFlow");
    if(!list) { printf("List not loaded\n"); return; }

    hist = (TH1D*) list->FindObject("fhRefsMult");
    hist->Scale(1/hist->GetEntries());
    list_fhRefsMult->Add(hist);
    hist = (TH1D*) list->FindObject("fhRefsPt");
    hist->Scale(1/hist->GetEntries());
    list_fhRefsPt->Add(hist);
    hist = (TH1D*) list->FindObject("fhRefsEta");
    hist->Scale(1/hist->GetEntries());
    list_fhRefsEta->Add(hist);
    hist = (TH1D*) list->FindObject("fhRefsPhi");
    hist->Scale(1/hist->GetEntries());
    list_fhRefsPhi->Add(hist);

    // V0s
    list = (TList*) gDirectory->Get("QA_V0s_UniFlow");
    if(!list) { printf("List not loaded\n"); return; }

    hist = (TH1D*) list->FindObject("fhQAV0sInvMassK0s_After");
    hist->Scale(1/hist->GetEntries());
    if(hist) list_fhQAV0sInvMassK0s->Add(hist);
    hist = (TH1D*) list->FindObject("fhQAV0sInvMassLambda_After");
    hist->Scale(1/hist->GetEntries());
    if(hist) list_fhQAV0sInvMassLambda->Add(hist);

    // Phi
    list = (TList*) gDirectory->Get("QA_Phi_UniFlow");
    if(!list) { printf("List not loaded\n"); return; }

    hist = (TH1D*) list->FindObject("fhPhiInvMass");
    hist->Scale(1/hist->GetEntries());
    if(hist) list_fhPhiInvMass->Add(hist);

    // FLOW histograms
    file = new TFile(fsFileNameFlow[iFile].Data(),"READ");
    if(!file) { printf("Input file '%s' not loaded\n",fsFileNameFlow[iFile].Data()); return; }

    // Refs
    file->ls();
    hist = (TH1D*) file->Get("hFlow2_Refs_harm2_gap08");
    list_hFlow2_Refs_harm2_gap08->Add(hist);

    // Charged
    hist = dynamic_cast<TH1D*>(file->Get("hFlow2_Charged_harm2_gap08_cent0;1"));
    if(hist) { hist->SetName(Form("%s_pos",hist->GetName())); list_hFlow2_Charged_harm2_gap08_cent0->Add(hist); }
    hist = dynamic_cast<TH1D*>(file->Get("hFlow2_Charged_harm2_gap08_cent0;2"));
    if(hist) { hist->SetName(Form("%s_neg",hist->GetName())); list_hFlow2_Charged_harm2_gap08_cent0->Add(hist); }

    hist = dynamic_cast<TH1D*>(file->Get("hFlow2_Charged_harm2_gap08_cent1;1"));
    if(hist) { hist->SetName(Form("%s_neg",hist->GetName())); list_hFlow2_Charged_harm2_gap08_cent1->Add(hist); }
    hist = dynamic_cast<TH1D*>(file->Get("hFlow2_Charged_harm2_gap08_cent1;2"));
    if(hist) { hist->SetName(Form("%s_neg",hist->GetName())); list_hFlow2_Charged_harm2_gap08_cent1->Add(hist); }

    hist = dynamic_cast<TH1D*>(file->Get("hFlow2_Charged_harm2_gap08_cent2;1"));
    if(hist) { hist->SetName(Form("%s_neg",hist->GetName())); list_hFlow2_Charged_harm2_gap08_cent2->Add(hist); }
    hist = dynamic_cast<TH1D*>(file->Get("hFlow2_Charged_harm2_gap08_cent2;2"));
    if(hist) { hist->SetName(Form("%s_neg",hist->GetName())); list_hFlow2_Charged_harm2_gap08_cent2->Add(hist); }

    hist = dynamic_cast<TH1D*>(file->Get("hFlow2_Charged_harm2_gap08_cent3;1"));
    if(hist) { hist->SetName(Form("%s_neg",hist->GetName())); list_hFlow2_Charged_harm2_gap08_cent3->Add(hist); }
    hist = dynamic_cast<TH1D*>(file->Get("hFlow2_Charged_harm2_gap08_cent3;2"));
    if(hist) { hist->SetName(Form("%s_neg",hist->GetName())); list_hFlow2_Charged_harm2_gap08_cent3->Add(hist); }

    hist = dynamic_cast<TH1D*>(file->Get("hFlow2_Charged_harm2_gap08_cent4;1"));
    if(hist) { hist->SetName(Form("%s_neg",hist->GetName())); list_hFlow2_Charged_harm2_gap08_cent4->Add(hist); }
    hist = dynamic_cast<TH1D*>(file->Get("hFlow2_Charged_harm2_gap08_cent4;2"));
    if(hist) { hist->SetName(Form("%s_neg",hist->GetName())); list_hFlow2_Charged_harm2_gap08_cent4->Add(hist); }
  }
  TCanvas* can = 0x0;

  can = PlotList(list_fhEventCentrality,1);
  can->SaveAs(Form("%s/fhEventCentrality.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotList(list_fhQAV0sInvMassK0s,1,"hist p");
  can->SaveAs(Form("%s/fhQAV0sInvMassK0s.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotList(list_fhQAV0sInvMassLambda,1,"hist p");
  can->SaveAs(Form("%s/fhQAV0sInvMassLambda.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotList(list_fhPhiInvMass,1, "hist p");
  can->SaveAs(Form("%s/fhPhiInvMass.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotList(list_fhRefsMult,1);
  can->SaveAs(Form("%s/fhRefsMult.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotList(list_fhRefsPt,1);
  can->SaveAs(Form("%s/fhRefsPt.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotList(list_fhRefsPhi,1);
  can->SaveAs(Form("%s/fhRefsPhi.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotList(list_fhRefsEta,1);
  can->SaveAs(Form("%s/fhRefsEta.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));

  can = PlotList(list_hFlow2_Refs_harm2_gap08,1,"hist p e");
  can->SaveAs(Form("%s/hFlow2_Refs_harm2_gap08.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotListPosNeg(list_hFlow2_Charged_harm2_gap08_cent0,1,"hist p e");
  can->SaveAs(Form("%s/hFlow2_Charged_harm2_gap08_cent0.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotListPosNeg(list_hFlow2_Charged_harm2_gap08_cent1,1,"hist p e");
  can->SaveAs(Form("%s/hFlow2_Charged_harm2_gap08_cent1.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotListPosNeg(list_hFlow2_Charged_harm2_gap08_cent2,1,"hist p e");
  can->SaveAs(Form("%s/hFlow2_Charged_harm2_gap08_cent2.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotListPosNeg(list_hFlow2_Charged_harm2_gap08_cent3,1,"hist p e");
  can->SaveAs(Form("%s/hFlow2_Charged_harm2_gap08_cent3.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  can = PlotListPosNeg(list_hFlow2_Charged_harm2_gap08_cent4,1,"hist p e");
  can->SaveAs(Form("%s/hFlow2_Charged_harm2_gap08_cent4.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));

  // TCanvas* c_hEventCentrality = PlotList(list_hEventCentrality,1);
  // c_hEventCentrality->SaveAs(Form("%s/hEventCentrality.%s",fsOutputPath.Data(),fsOutputFileFormat.Data()));
  // TCanvas* c_hEventCentrality = PlotList(list_hEventCentrality,1);
  return;
}
// _____________________________________________________________________________
TCanvas* PlotList(const TList* list, const Short_t iPlotRatio, const char* sDrawOpt)
{
  printf(" ------- PlotList -------\n");
  if(!list) { printf("Input list does not exists!\n"); return 0x0; }

  const Short_t iNumEntries = list->GetEntries();

  printf("  Entries: %d\n",iNumEntries);
  if(iNumEntries == 0) { printf("Input list is empty!\n"); return 0x0; }

  TH1* histNom = list->At(0); // placeholder for ratio nominator (common)
  TH1* histRatio = 0x0; // placeholder for ratio
  TH1* hist = 0x0; // placeholder for temporary histo

  TLegend* legend = new TLegend(0.6,0.6,0.8,0.8);

  TCanvas* canvas = new TCanvas(histNom->GetName(),histNom->GetName(),800,800);
  if(iPlotRatio)
  {
    canvas->Divide(1,2);
    // TH1* cFrameRatio = canvas->cd(2)->DrawFrame(0,0,100,2);
  }

  // finding the maximum and minimum y-value
  Double_t dMax = 0;
  Double_t dMin = 999999;
  for(Short_t iFile(0); iFile < iNumEntries; iFile++)
  {
    hist = (list->At(iFile)->Class()) list->At(iFile);
    if(!hist) { printf("File %d: input histo not loaded properly\n",iFile); return 0x0; }
    if(dMax < hist->GetMaximum()) dMax = hist->GetMaximum();
    if(dMin > hist->GetMinimum()) dMin = hist->GetMinimum();
  }

  for(Short_t iFile(0); iFile < fiNumFiles; iFile++)
  {
    // printf("Class: '%s'\n", list->At(0)->ClassName());
    hist = (list->At(iFile)->Class()) list->At(iFile);
    if(!hist) { printf("File %d: input histo not loaded properly\n",iFile); return 0x0; }

    hist->SetMaximum(dMax*1.3);
    hist->SetMinimum(dMin*0.7);
    hist->SetStats(0);
    hist->SetMarkerStyle(fMarkers[iFile]);
    hist->SetMarkerColor(fColors[iFile]);
    hist->SetLineColor(fColors[iFile]);

    legend->AddEntry(hist,fsTriggerName[iFile].Data(),"pl");

    canvas->cd(1);
    (iFile == 0) ? hist->Draw(sDrawOpt) : hist->Draw(Form("%s same",sDrawOpt));

    if(iPlotRatio && iFile > 0)
    {
      canvas->cd(2);
      histRatio = (list->At(iFile)->Class()) histNom->Clone();
      histRatio->SetTitle("");
      histRatio->Divide(hist);
      histRatio->SetMaximum(2.);
      histRatio->SetMinimum(0.);
      histRatio->SetStats(0);
      histRatio->SetMarkerStyle(fMarkers[iFile]);
      histRatio->SetMarkerColor(fColors[iFile]);
      histRatio->SetLineColor(fColors[iFile]);
      histRatio->Draw("same");
    }
  }

  canvas->cd(1);
  legend->Draw();

  if(iPlotRatio)
  {
    canvas->cd(2);
    TLine* unity = new TLine();
    unity->SetLineColor(fColors[0]);
    unity->DrawLine(histRatio->GetXaxis()->GetXmin(),1.,histRatio->GetXaxis()->GetXmax() ,1.);
  }

  return canvas;
}
// _____________________________________________________________________________
TCanvas* PlotListPosNeg(const TList* list, const Short_t iPlotRatio, const char* sDrawOpt)
{
  printf(" ------- PlotListPosNeg -------\n");
  if(!list) { printf("Input list does not exists!\n"); return 0x0; }

  const Short_t iNumEntries = list->GetEntries();

  printf("  Entries: %d\n",iNumEntries);
  if(iNumEntries == 0) { printf("Input list is empty!\n"); return 0x0; }

  // list->ls();

  TH1* histNom = list->At(0); // placeholder for ratio nominator (common)
  TH1* histNom2 = list->At(1); // placeholder for ratio nominator (common)
  TH1* histRatio = 0x0; // placeholder for ratio
  TH1* hist = 0x0; // placeholder for temporary histo
  TH1* hist2 = 0x0; // placeholder for temporary histo

  const Double_t dRatioYmax = 1.5;
  const Double_t dRatioYmin = 0.5;

  // TLegend* legend = new TLegend(0.12,0.5,0.4,0.89);
  TLegend* legend = new TLegend(0.76,0.3,0.99,0.89);
  legend->SetBorderSize(0);
  // legend->SetFillColorAlpha(0,0);

  TLegend* legendRatio = new TLegend(0.76,0.25,0.99,0.6);
  legendRatio->SetBorderSize(0);
  // legendRatio->SetFillColorAlpha(0,0);

  TLegend* legendRatioPosNeg = new TLegend(0.76,0.25,0.99,0.6);
  legendRatioPosNeg->SetBorderSize(0);
  // legendRatioPosNeg->SetFillColorAlpha(0,0);


  TCanvas* canvas = new TCanvas(histNom->GetName(),histNom->GetName(),1200,1000);
  if(iPlotRatio)
  {
    canvas->Divide(1,3);
    // TH1* cFrameRatio = canvas->cd(2)->DrawFrame(0,0,100,2);

    TLine* unity = new TLine();
    unity->SetLineColor(fColors[0]);

  }
  canvas->cd(1);
  gPad->SetRightMargin(0.25);
  canvas->cd(2);
  gPad->SetRightMargin(0.25);
  canvas->cd(3);
  gPad->SetRightMargin(0.25);


  // finding the maximum and minimum y-value
  Double_t dMax = 0;
  Double_t dMin = 999999;
  for(Short_t iFile(0); iFile < iNumEntries; iFile++)
  {
    hist = (list->At(iFile)->Class()) list->At(iFile);
    if(!hist) { printf("File %d: input histo not loaded properly\n",iFile); return 0x0; }
    if(dMax < hist->GetMaximum()) dMax = hist->GetMaximum();
    if(dMin > hist->GetMinimum()) dMin = hist->GetMinimum();
  }

  for(Short_t iFile(0); iFile < fiNumFiles; iFile++)
  {
    // printf("Class: '%s'\n", list->At(0)->ClassName());
    hist = (list->At(iFile)->Class()) list->At(2*iFile);
    hist2 = (list->At(iFile)->Class()) list->At(2*iFile+1);
    if(!hist) { printf("File %d: input histo not loaded properly\n",iFile); return 0x0; }
    if(!hist2) { printf("File %d: input histo not loaded properly\n",iFile); return 0x0; }

    hist->SetMaximum(dMax*1.3);
    hist->SetMinimum(dMin*0.7);
    hist->SetStats(0);
    hist->SetMarkerStyle(fMarkers[iFile]);
    hist->SetMarkerColor(fColors[iFile]);
    hist->SetLineColor(fColors[iFile]);

    hist2->SetMaximum(dMax*1.3);
    hist2->SetMinimum(dMin*0.7);
    hist2->SetStats(0);
    hist2->SetMarkerStyle(fMarkers[iFile+3]);
    hist2->SetMarkerColor(fColors[iFile]);
    hist2->SetLineColor(fColors[iFile]);

    legend->AddEntry(hist,Form("%s (POI pos)",fsTriggerName[iFile].Data()),"pl");
    legend->AddEntry(hist2,Form("%s (POI neg)",fsTriggerName[iFile].Data()),"pl");

    canvas->cd(1);
    (iFile == 0) ? hist->Draw(sDrawOpt) : hist->Draw(Form("%s same",sDrawOpt));
    hist2->Draw(Form("%s same",sDrawOpt));

    if(iPlotRatio)
    {
      if(iFile > 0)
      {
        canvas->cd(2);
        histRatio = (list->At(iFile)->Class()) histNom->Clone();
        histRatio->SetTitle("Ratio: FAST / CENT (w/wo SDD)");
        histRatio->Divide(hist);
        histRatio->SetMaximum(dRatioYmax);
        histRatio->SetMinimum(dRatioYmin);
        histRatio->SetStats(0);
        histRatio->SetMarkerStyle(fMarkers[iFile]);
        histRatio->SetMarkerColor(fColors[iFile]);
        histRatio->SetLineColor(fColors[iFile]);
        histRatio->Draw("same");
        legendRatio->AddEntry(histRatio,Form("POI pos: FAST / %s",fsTriggerName[iFile].Data()),"pl");

        histRatio = (list->At(iFile)->Class()) histNom2->Clone();
        histRatio->Divide(hist);
        histRatio->SetMaximum(dRatioYmax);
        histRatio->SetMinimum(dRatioYmin);
        histRatio->SetStats(0);
        histRatio->SetMarkerStyle(fMarkers[iFile+3]);
        histRatio->SetMarkerColor(fColors[iFile]);
        histRatio->SetLineColor(fColors[iFile]);
        histRatio->Draw("same");
        legendRatio->AddEntry(histRatio,Form("POI neg: FAST / %s",fsTriggerName[iFile].Data()),"pl");

        if(iFile == 1) unity->DrawLine(histRatio->GetXaxis()->GetXmin(),1.,histRatio->GetXaxis()->GetXmax() ,1.);
      }

      canvas->cd(3);
      histRatio = (list->At(iFile)->Class()) hist->Clone();
      histRatio->SetTitle("Ration: POIs positive eta / negative eta");
      histRatio->Divide(hist2);
      histRatio->SetMaximum(dRatioYmax);
      histRatio->SetMinimum(dRatioYmin);
      histRatio->SetStats(0);
      histRatio->SetMarkerStyle(fMarkers[iFile]);
      histRatio->SetMarkerColor(fColors[iFile]);
      histRatio->SetLineColor(fColors[iFile]);
      histRatio->Draw("same");
      legendRatioPosNeg->AddEntry(histRatio,Form("%s: POI pos / POI neg",fsTriggerName[iFile].Data()),"pl");
      if(iFile == 1) unity->DrawLine(histRatio->GetXaxis()->GetXmin(),1.,histRatio->GetXaxis()->GetXmax() ,1.);

    }
  }

  canvas->cd(1);
  legend->Draw();
  canvas->cd(2);
  legendRatio->Draw();
  canvas->cd(3);
  legendRatioPosNeg->Draw();

  return canvas;
}
// _____________________________________________________________________________
void PlotComparison(TString histName, TString histListName, const char* histClass)
{
  printf("====== PlotComparison: histName '%s' | histListName '%s' ======\n",histName.Data(),histListName.Data());
  if (histName.EqualTo("")) { printf("Hist name is empty.\n"); return; }
  if (histListName.EqualTo("")) { printf("Hist name is empty.\n"); return; }
  if (histListName.EqualTo("")) { printf("Hist class is empty.\n"); return; }

  TLegend* legend = new TLegend(0.4,0.4,0.8,0.8);

  TCanvas* canvas = new TCanvas(histName.Data(),histName.Data(),800,800);
  canvas->cd();

  TFile* file = 0x0;
  TList* list = 0x0;
  TH1* hist = 0x0;
  for(Short_t iFile(0); iFile < fiNumFiles; iFile++)
  {
    file = new TFile(fsFileName[iFile].Data(),"READ");
    if(!file) { printf("No input file: '%s'\n",fsFileName[iFile].Data()); return; }
    file->cd("UniFlow");
    // file->ls();

    list = (TList*) gDirectory->Get(histListName.Data());
    if(!list) {printf("File '%s': Input hist list '%s' not found! Check its name.\n",fsFileName[iFile],histListName.Data()); return; }

    switch(histClass)
    {
      case "TH1D":
        hist = (TH1D*) list->FindObject(histName.Data());
      break;

      case "TH2D":
        hist = (TH2D*) list->FindObject(histName.Data());
      break;

      case "TH3D":
        hist = (TH3D*) list->FindObject(histName.Data());
      break;

      default:
        printf("hist class '%s' not supported.\n",histClass);
        return;
    }

    if(!hist) { printf("File %s: Hist '%s' not loaded properly.",fsFileName[iFile].Data(), histName.Data()); return; }

    hist->SetMarkerStyle(fMarkers[iFile]);
    hist->SetMarkerColor(fColors[iFile]);
    hist->SetLineColor(fColors[iFile]);

    legend->AddEntry(hist,fsTriggerName[iFile].Data(),"pl");

    (iFile == 0) ? hist->Draw() : hist->Draw("same");

  }

  legend->Draw();

  return;
}
