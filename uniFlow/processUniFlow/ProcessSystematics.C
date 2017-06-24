/* ProcessSystematics
 *
 * Macro for processing results of systematical test.
 * Producing the ratios and estimating the overall uncertainty.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

#include "TROOT.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"

const char* sInputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2";
const char* sOutputPath = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/comparison";

const char* sBaseFile = "/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTest/UniFlow_NUA_cor.root";

// const char* sCommonTag = "tpcCls"; const Short_t iNumFiles = 1; const char* sTaskTag[iNumFiles] = {"uniflow_tpcCls80"};
const char* sCommonTag = "tpcCls"; const Short_t iNumFiles = 2; const char* sTaskTag[iNumFiles] = {"uniflow_tpcCls80","uniflow_tpcCls90"};
// const char* sTaskTag[iNumFiles] = {"NUA_cor_cent","noNUA_cent","vtx8","vtx9","tpcCls80","tpcCls90","CPA","filterBit","bgFunc","fitRange"};
// const char* sTaskTag[iNumFiles] = {"NUA_cor_cent","noNUA_cent","vtx8","vtx9","tpcCls80","tpcCls90","2sigma","bayes","filterBit","bgFunc","fitRange"};
// const char* sTaskTag[iNumFiles] = {"NUA_cor_cent","noNUA_cent","vtx8","vtx9","tpcCls80","tpcCls90","filterBit"};

// const char* sTaskTag[iNumFiles] = {"fb96","fb768"};

//const Short_t iNumSpecies = 4; const char* sSpecies[iNumSpecies] = {"hFlow2_Charged_harm2_gap08_cent","hFlow2_Pion_harm2_gap08_cent","hFlow2_Kaon_harm2_gap08_cent","hFlow2_Proton_harm2_gap08_cent"};
// const Short_t iNumSpecies = 3; const char* sSpecies[iNumSpecies] = {"hFlow2_K0s_harm2_gap08_mult","hFlow2_Lambda_harm2_gap08_mult","hFlow2_Phi_harm2_gap08_mult"};
const Short_t iNumSpecies = 2; const char* sSpecies[iNumSpecies] = {"hFlow2_K0s_harm2_gap08_mult","hFlow2_Lambda_harm2_gap08_mult"}; Double_t dFittRangeLow[iNumSpecies] = {0.5, 0.5}; Double_t dFittRangeHigh[iNumSpecies] = {4.,4.};
// const Short_t iNumSpecies = 1; const char* sSpecies[iNumSpecies] = {"hFlow2_Phi_harm2_gap08_mult"};



// Bool_t plotTogether = kTRUE;
Bool_t plotTogether = kFALSE;

const Short_t iNumCent = 4;
const char* sCentTag[iNumCent] = {"0-20","20-40","40-60","60-100"};

const Int_t fillColors[] = {kGray+1,  kRed-10, kBlue-9, kGreen-8, kMagenta-9, kOrange-9,kCyan-8,kYellow-7};
const Int_t colors[]     = {kGreen+3, kBlue+1,kRed,kOrange-1,kCyan+2,kYellow+2,kRed+2,kMagenta-3,kBlack};
const Int_t markers[]    = {kFullCircle, kFullSquare,kFullDiamond,kFullCross,kOpenCircle,kOpenSquare,kOpenDiamond,kOpenCross,kFullDiamond,kFullStar,kOpenStar};



// TLegend* legend;

// void ProcessList_1by1(TList* list = 0x0);
// void ProcessList(TList* list = 0x0);
TH1D* MakeRatio(TH1D* hBase = 0x0, TH1D* hSyst = 0x0, Bool_t bCorelated = kFALSE);
TH1D* DoBarlowTest(TH1D* hBase = 0x0, TH1D* hSyst = 0x0, Bool_t bCorelated = kFALSE);

void ProcessSystematics()
{
  // // legend = new TLegend(0.62, 0.18, 0.8, 0.46);

  gStyle->SetOptFit(1111);
  // gStyle->SetOptStat(1000);

  TList* list = new TList();

  TCanvas* canSuper = new TCanvas("canSuper","canSuper",800,800);
  canSuper->Divide(3,4);
  TCanvas* canRatio = new TCanvas("canRatio","canRatio",800,800);
  TCanvas* canBarlow = new TCanvas("canBarlow","canBarlow",800,800);
  TCanvas* canOver = new TCanvas("canOver","canOver",1200,400);
  canOver->Divide(3,1);

  TH1* frame_ratio = 0x0;
  TH1* frame_barlow = 0x0;

  TFile* fInputFile = 0x0;
  TH1D* hBase = 0x0;
  TH1D* hTemp = 0x0;
  TH1D* hRatio = 0x0;
  TH1D* hBarlow = 0x0;

  TF1* fitRatio = 0x0;

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.06);
  TLatex* latexRed = new TLatex();
  latexRed->SetTextColor(kRed);
  latexRed->SetNDC();
  latexRed->SetTextSize(0.06);

  TLatex* latexTog = new TLatex();
  latexTog->SetNDC();
  latexTog->SetTextSize(0.06);

  TFile* fInputFileBase = new TFile(sBaseFile,"READ");
  if(!fInputFileBase) { printf("InputFileBase not found!\n"); return; }

  TLegend* legend = new TLegend(0.45,0.14,0.82,0.55);
  TLegend* legendSolo = new TLegend(0.45,0.25,0.82,0.45);
  legendSolo->SetBorderSize(0);
  legendSolo->SetFillColorAlpha(0,0);
  legend->SetBorderSize(0);
  legend->SetFillColorAlpha(0,0);


  if(!plotTogether)
  {
    for(Short_t iFile(0); iFile < iNumFiles; iFile++)
    {
      fInputFile = new TFile(Form("%s/%s/%s/UniFlow_%s.root",sInputPath,sCommonTag,sTaskTag[iFile],sTaskTag[iFile]),"READ");
      if(!fInputFile) { printf("InputFile (%d) not found!\n",iFile); return; }

      // loop over species (histos of interests)
      for(Short_t iSpecies(0); iSpecies < iNumSpecies; iSpecies++)
      {
        // loop over centralities (histos of interests)
        for(Short_t iCent(0); iCent < iNumCent; iCent++)
        {

          hBase = (TH1D*) fInputFileBase->Get(Form("%s%d",sSpecies[iSpecies],iCent));
          if(!hBase) { printf("Histos Base not found: file %d | species %d\n",iFile,iSpecies); fInputFileBase->ls(); return; }
          hBase->SetLineColor(colors[0]);
          hBase->SetMarkerColor(colors[0]);
          hBase->SetMarkerStyle(markers[0]);
          hBase->SetTitle("; #it{p}_{T} (GeV#it{c}); v_{2}");

          hTemp = (TH1D*) fInputFile->Get(Form("%s%d",sSpecies[iSpecies],iCent));
          if(!hTemp) { printf("Histos not found: file %d | species %d\n",iFile,iSpecies); fInputFile->ls(); return; }
          hTemp->SetLineColor(colors[1]);
          hTemp->SetMarkerColor(colors[1]);
          hTemp->SetMarkerStyle(markers[1]);

          if(iCent == 0 && iSpecies == 0) legendSolo->AddEntry(hTemp,sTaskTag[iFile],"pel");


          hRatio = MakeRatio(hBase,hTemp,1);
          hRatio->SetTitle("Ratio systematic/baseline; #it{p}_{T} (GeV/#it{c}); syst / baseline;");
          hRatio->SetName(Form("default/%s",sTaskTag[iFile]));

          hBarlow = DoBarlowTest(hBase,hTemp,1);
          hBarlow->SetTitle("Barlow Test; #it{p}_{T} (GeV/#it{c}); Barlow #it{B}");
          hBarlow->SetMaximum(5);
          hBarlow->SetMinimum(0);

          canRatio->cd();
          gPad->SetGrid(0,1);
          // frame_ratio = gPad->DrawFrame(hBase->GetXaxis()->GetXmin(),0.,hBase->GetXaxis()->GetXmax(),2.);
          hRatio->Draw("hist p e ");
          // fitRatio = new TF1("fitRatio","[0]",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
          fitRatio = new TF1("fitRatio","[0]",dFittRangeLow[iSpecies],dFittRangeHigh[iSpecies]);
          hRatio->Fit("fitRatio","R");
          fitRatio->Draw("same");
          latex->DrawLatex(0.16,0.77,Form("%s%%",sCentTag[iCent]));
          latexRed->DrawLatex(0.45,0.77,Form("%.3g#pm%.3g",fitRatio->GetParameter(0),fitRatio->GetParError(0)));
          latex->DrawLatex(0.16,0.84,Form("%s",sSpecies[iSpecies]));
          // canRatio->SaveAs(Form("%s/%s/Ratio_%s_%s_cent%d.pdf",sOutputPath,sCommonTag,sSpecies[iSpecies],sTaskTag[iFile],iCent));

          canBarlow->cd();
          gPad->SetGrid(0,1);
          // frame_barlow = gPad->DrawFrame(hBase->GetXaxis()->GetXmin(),0.,hBase->GetXaxis()->GetXmax(),5.);
          hBarlow->Draw("hist p e");
          // canBarlow->SaveAs(Form("%s/%s/Barlow_%s_%s_cent%d.pdf",sOutputPath,sCommonTag,sSpecies[iSpecies],sTaskTag[iFile],iCent));

          canOver->cd(1);
          hBase->Draw("hist p e");
          hTemp->Draw("hist p e same");
          canOver->cd(2);
          gPad->SetGrid(0,1);
          // frame_ratio = gPad->DrawFrame(hBase->GetXaxis()->GetXmin(),0.,hBase->GetXaxis()->GetXmax(),2.);
          hRatio->Draw("hist p e");
          fitRatio->Draw("same");
          latex->DrawLatex(0.13,0.74,Form("%s%%",sCentTag[iCent]));
          latex->DrawLatex(0.13,0.82,Form("%s",sSpecies[iSpecies]));
          latexRed->DrawLatex(0.45,0.74,Form("%.3g#pm%.3g",fitRatio->GetParameter(0),fitRatio->GetParError(0)));
          canOver->cd(3);
          // frame_barlow = gPad->DrawFrame(hBase->GetXaxis()->GetXmin(),0.,hBase->GetXaxis()->GetXmax(),5.);
          gPad->SetGrid(0,1);
          hBarlow->Draw("hist p e");
          // canOver->SaveAs(Form("%s/%s/Over_%s_%s_cent%d.pdf",sOutputPath,sCommonTag,sSpecies[iSpecies],sTaskTag[iFile],iCent));

          canSuper->cd(iCent*3+1);
          // gPad->SetGrid(0,1);
          hBase->Draw("hist p e");
          hTemp->Draw("hist p e same");
          // if(iCent == 0)
          legendSolo->Draw();



          canSuper->cd(iCent*3+2);
          gPad->SetGrid(0,1);
          hRatio->Draw("hist p e");
          fitRatio->Draw("same");
          latex->DrawLatex(0.16,0.77,Form("%s%%",sCentTag[iCent]));
          latexRed->DrawLatex(0.45,0.77,Form("%.3g#pm%.3g",fitRatio->GetParameter(0),fitRatio->GetParError(0)));
          latex->DrawLatex(0.16,0.84,Form("%s",sSpecies[iSpecies]));

          canSuper->cd(iCent*3+3);
          gPad->SetGrid(0,1);
          hBarlow->Draw("hist p e");
        }
        canSuper->SaveAs(Form("%s/%s/0_%s_%s.pdf",sOutputPath,sCommonTag,sSpecies[iSpecies],sTaskTag[iFile]));
      }
    }
  }
  else
  {
    // loop over species (histos of interests)
    for(Short_t iSpecies(0); iSpecies < iNumSpecies; iSpecies++)
    {
      for(Short_t iFile(0); iFile < iNumFiles; iFile++)
      {
        // canSuper->Clear("D");
        fInputFile = new TFile(Form("%s/%s/%s/UniFlow_%s.root",sInputPath,sCommonTag,sTaskTag[iFile],sTaskTag[iFile]),"READ");
        if(!fInputFile) { printf("InputFile (%d) not found!\n",iFile); return; }

        // loop over centralities (histos of interests)
        for(Short_t iCent(0); iCent < iNumCent; iCent++)
        {
          hBase = (TH1D*) fInputFileBase->Get(Form("%s%d",sSpecies[iSpecies],iCent));
          if(!hBase) { printf("Histos Base not found: file %d | species %d\n",iFile,iSpecies); fInputFileBase->ls(); return; }
          hBase->SetLineColor(colors[0]);
          hBase->SetMarkerColor(colors[0]);
          hBase->SetMarkerStyle(markers[0]);
          hBase->SetTitle("; #it{p}_{T} (GeV#it{c}); v_{2}");

          hTemp = (TH1D*) fInputFile->Get(Form("%s%d",sSpecies[iSpecies],iCent));
          if(!hTemp) { printf("Histos not found: file %d | species %d\n",iFile,iSpecies); fInputFile->ls(); return; }
          hTemp->SetLineColor(colors[iFile+1]);
          hTemp->SetMarkerColor(colors[iFile+1]);
          hTemp->SetMarkerStyle(markers[iFile+1]);

          if(iCent == 0 && iSpecies == 0) legend->AddEntry(hTemp,sTaskTag[iFile],"pel");


          hRatio = MakeRatio(hBase,hTemp,1);
          hRatio->SetTitle("Ratio systematic/baseline; #it{p}_{T} (GeV/#it{c}); syst / baseline;");
          hRatio->SetName(Form("default/%s",sTaskTag[iFile]));

          hBarlow = DoBarlowTest(hBase,hTemp,1);
          hBarlow->SetTitle("Barlow Test; #it{p}_{T} (GeV/#it{c}); Barlow #it{B}");
          hBarlow->SetMaximum(5);
          hBarlow->SetMinimum(0);

          canRatio->cd();
          gPad->SetGrid(0,1);
          // frame_ratio = gPad->DrawFrame(hBase->GetXaxis()->GetXmin(),0.,hBase->GetXaxis()->GetXmax(),2.);
          hRatio->Draw("hist p e ");
          // fitRatio = new TF1("fitRatio","[0]",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
          fitRatio = new TF1("fitRatio","[0]",dFittRangeLow[iSpecies],dFittRangeHigh[iSpecies]);
          fitRatio->SetLineColor(colors[iFile+1]);
          hRatio->Fit("fitRatio","R");
          fitRatio->Draw("same");
          latex->DrawLatex(0.16,0.77,Form("%s%%",sCentTag[iCent]));
          latexRed->DrawLatex(0.45,0.77,Form("%.3g#pm%.3g",fitRatio->GetParameter(0),fitRatio->GetParError(0)));
          latex->DrawLatex(0.16,0.84,Form("%s",sSpecies[iSpecies]));
          // canRatio->SaveAs(Form("%s/%s/Ratio_%s_%s_cent%d.pdf",sOutputPath,sCommonTag,sSpecies[iSpecies],sTaskTag[iFile],iCent));

          canBarlow->cd();
          gPad->SetGrid(0,1);
          // frame_barlow = gPad->DrawFrame(hBase->GetXaxis()->GetXmin(),0.,hBase->GetXaxis()->GetXmax(),5.);
          hBarlow->Draw("hist p e");
          // canBarlow->SaveAs(Form("%s/%s/Barlow_%s_%s_cent%d.pdf",sOutputPath,sCommonTag,sSpecies[iSpecies],sTaskTag[iFile],iCent));

          canOver->cd(1);
          hBase->Draw("hist p e");
          hTemp->Draw("hist p e same");
          canOver->cd(2);
          gPad->SetGrid(0,1);
          // frame_ratio = gPad->DrawFrame(hBase->GetXaxis()->GetXmin(),0.,hBase->GetXaxis()->GetXmax(),2.);
          hRatio->Draw("hist p e");
          fitRatio->Draw("same");
          latex->DrawLatex(0.13,0.74,Form("%s%%",sCentTag[iCent]));
          latex->DrawLatex(0.13,0.82,Form("%s",sSpecies[iSpecies]));
          latexRed->DrawLatex(0.45,0.74,Form("%.3g#pm%.3g",fitRatio->GetParameter(0),fitRatio->GetParError(0)));
          canOver->cd(3);
          // frame_barlow = gPad->DrawFrame(hBase->GetXaxis()->GetXmin(),0.,hBase->GetXaxis()->GetXmax(),5.);
          gPad->SetGrid(0,1);
          hBarlow->Draw("hist p e");
          // canOver->SaveAs(Form("%s/%s/Over_%s_%s_cent%d.pdf",sOutputPath,sCommonTag,sSpecies[iSpecies],sTaskTag[iFile],iCent));

          canSuper->cd(iCent*3+1);
          // gPad->SetGrid(0,1);
          if(iFile == 0) hBase->Draw("hist p e ");
          else hBase->Draw("hist p e same");
          hTemp->Draw("hist p e same");
          if(iFile == 0)legend->Draw();


          canSuper->cd(iCent*3+2);
          gPad->SetGrid(0,1);
          if(iFile == 0) hRatio->Draw("hist p e ");
          else hRatio->Draw("same hist p e");
          fitRatio->Draw("same");
          latex->DrawLatex(0.16,0.77,Form("%s%%",sCentTag[iCent]));
          latexTog->SetTextColor(colors[iFile+1]);
          latexTog->DrawLatex(0.45,0.17+iFile*0.1,Form("%.3g#pm%.3g",fitRatio->GetParameter(0),fitRatio->GetParError(0)));
          latex->DrawLatex(0.16,0.84,Form("%s",sSpecies[iSpecies]));

          canSuper->cd(iCent*3+3);
          gPad->SetGrid(0,1);
          if(iFile == 0) hBarlow->Draw(" hist p e");
          else hBarlow->Draw(" hist p e same");
        }
        canSuper->SaveAs(Form("%s/%s/0_%s_%s.pdf",sOutputPath,sCommonTag,sCommonTag,sSpecies[iSpecies]));
      }
    }

  }


  // // loop over species (histos of interests)
  // for(Short_t iSpecies(0); iSpecies < iNumSpecies; iSpecies++)
  // {
  //   // loop over centralities (histos of interests)
  //   for(Short_t iCent(0); iCent < iNumCent; iCent++)
  //   {
  //     // loop over files
  //     for(Short_t iFile(0); iFile < iNumFiles; iFile++)
  //     {
  //       // fInputFile = new TFile(Form("%s/%s/UniFlow_%s.root",sInputPath,sTaskTag[iFile],sTaskTag[iFile]),"READ");
  //       fInputFile = new TFile(Form("%s/UniFlow_%s.root",sInputPath,sTaskTag[iFile]),"READ");
  //       if(!fInputFile) { printf("InputFile (%d) not found!\n",iFile); return; }
  //
  //       hTemp = (TH1D*) fInputFile->Get(Form("%s%d",sSpecies[iSpecies],iCent));
  //       if(!hTemp) { printf("Histos not found: file %d | species %d\n",iFile,iSpecies); fInputFile->ls(); return; }
  //       list->Add(hTemp);
  //     }
  //
  //     ProcessList_1by1(list);
  //     // saving canvases
  //     // canSuper->SaveAs(Form("%s/%s_cent%d_sup.pdf",sOutputPath,sSpecies[iSpecies],iCent));
  //     canSuper->SaveAs(Form("%s/%s_cent%d_sup.png",sOutputPath,sSpecies[iSpecies],iCent));
  //     // canSuper->SaveAs(Form("%s/%s_cent%d_sup.eps",sOutputPath,sSpecies[iSpecies],iCent));
  //
  //     // canRatio->SaveAs(Form("%s/%s_cent%d_rat.pdf",sOutputPath,sSpecies[iSpecies],iCent));
  //     canRatio->SaveAs(Form("%s/%s_cent%d_rat.png",sOutputPath,sSpecies[iSpecies],iCent));
  //     // canRatio->SaveAs(Form("%s/%s_cent%d_rat.eps",sOutputPath,sSpecies[iSpecies],iCent));
  //
  //     // canBarlow->SaveAs(Form("%s/%s_cent%d_bar.pdf",sOutputPath,sSpecies[iSpecies],iCent));
  //     canBarlow->SaveAs(Form("%s/%s_cent%d_bar.png",sOutputPath,sSpecies[iSpecies],iCent));
  //     // canBarlow->SaveAs(Form("%s/%s_cent%d_bar.eps",sOutputPath,sSpecies[iSpecies],iCent));
  //
  //     canSuper->Clear();
  //     canRatio->Clear();
  //     canBarlow->Clear();
  //     list->Clear();
  //   }
  // }

  // printf("List entries: %d\n",list->GetEntries());

  return;
}
// //_____________________________________________________________________________
// void ProcessList(TList* list)
// {
//   if(!list) { printf("ProcessList::List not found\n"); return; }
//   const Short_t iEntries = list->GetEntries();
//   if(iEntries != iNumFiles) { printf("Different number of entries (%d) than expected (%d)\n",iEntries,iNumFiles); return; }
//
//   TH1D* hBase = (TH1D*) list->At(0);
//   if(!hBase) { printf("ProcessList::Baseline histo not found\n"); return; }
//
//   // legend with default
//   TLegend* legendFull = new TLegend(0.62, 0.18, 0.8, 0.46);
//   legendFull->SetBorderSize(0);
//   legendFull->SetFillColorAlpha(0,0);
//   legendFull->AddEntry(hBase,"default","pl");
//
//   // legend without default
//   TLegend* legend = new TLegend(0.62, 0.18, 0.8, 0.46);
//   legend->SetBorderSize(0);
//   legend->SetFillColorAlpha(0,0);
//
//   hBase->SetLineColor(colors[0]);
//   hBase->SetMarkerColor(colors[0]);
//   hBase->SetMarkerStyle(markers[0]);
//   hBase->SetFillColor(colors[0]);
//   hBase->SetMinimum(0.);
//   hBase->SetMaximum(1.);
//
//   canSuper->cd();
//   TH1* frame_super = canSuper->DrawFrame(0.,0.,8.,1.);
//   hBase->Draw("same");
//
//   canRatio->cd();
//   TH1* frame_ratio = canRatio->DrawFrame(0.,0.,8.,2.);
//   canBarlow->cd();
//   TH1* frame_barlow = canBarlow->DrawFrame(0.,0.,8.,10.);
//
//
//   TH1D* hTemp = 0x0;
//   TH1D* hRatio = 0x0;
//   // TH1D* hRatio2 = 0x0;
//   TH1D* hBarlow = 0x0;
//
//   for(Short_t iFile(1); iFile < iNumFiles; iFile++)
//   {
//     hTemp = (TH1D*) list->At(iFile);
//     if(!hTemp) { printf("ProcessList::Histo (%d) not found\n",iFile); return; }
//
//     hTemp->SetLineColor(colors[iFile]);
//     hTemp->SetMarkerColor(colors[iFile]);
//     hTemp->SetMarkerStyle(markers[iFile]);
//     hTemp->SetFillColor(colors[iFile]);
//
//     legendFull->AddEntry(hTemp,sTaskTag[iFile],"pl");
//     legend->AddEntry(hTemp,sTaskTag[iFile],"pl");
//
//     canSuper->cd();
//     hTemp->Draw("same");
//
//     hRatio = MakeRatio(hBase,hTemp,1);
//     // hRatio->SetLineColor(colors[iFile]);
//     // hRatio->SetMarkerColor(colors[iFile]);
//     // hRatio->SetMarkerStyle(markers[iFile]);
//     // hRatio->SetFillColor(colors[iFile]);
//
//     hBarlow = DoBarlowTest(hBase,hTemp,1);
//
//     canRatio->cd();
//     hRatio->Draw("hist p e1 same");
//
//     canBarlow->cd();
//     hBarlow->Draw("hist p same");
//   }
//
//   canSuper->cd();
//   legendFull->Draw();
//
//   // fitting the ratios
//   canRatio->cd();
//   legend->Draw();
//   TF1* fitConst = new TF1("fitConst","[0]",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
//   fitConst->SetLineColor(kRed);
//   hRatio->Fit("fitConst","R");
//   printf("Fit: %g ± %g (Chi2 %g)\n",fitConst->GetParameter(0),fitConst->GetParError(0),fitConst->GetChisquare());
//
//   TLine* lineUnity = new TLine();
//   lineUnity->SetLineColor(kBlack);
//   lineUnity->SetLineStyle(7);
//
//   fitConst->Draw("same");
//   lineUnity->DrawLine(hRatio->GetXaxis()->GetXmin(),1.,hRatio->GetXaxis()->GetXmax(),1.);
//
//   canBarlow->cd();
//   lineUnity->DrawLine(hBarlow->GetXaxis()->GetXmin(),1.,hBarlow->GetXaxis()->GetXmax(),1.);
//   legend->Draw();
//
//   return;
// }
// //_____________________________________________________________________________
// //_____________________________________________________________________________
// void ProcessList_1by1(TList* list)
// {
//   if(!list) { printf("ProcessList::List not found\n"); return; }
//   const Short_t iEntries = list->GetEntries();
//   if(iEntries != iNumFiles) { printf("Different number of entries (%d) than expected (%d)\n",iEntries,iNumFiles); return; }
//
//   TH1D* hBase = (TH1D*) list->At(0);
//   if(!hBase) { printf("ProcessList::Baseline histo not found\n"); return; }
//
//   // legend with default
//   TLegend* legendFull = new TLegend(0.62, 0.18, 0.8, 0.46);
//   legendFull->SetBorderSize(0);
//   legendFull->SetFillColorAlpha(0,0);
//   legendFull->AddEntry(hBase,"default","pl");
//
//   // legend without default
//   TLegend* legend = new TLegend(0.62, 0.18, 0.8, 0.46);
//   legend->SetBorderSize(0);
//   legend->SetFillColorAlpha(0,0);
//
//   hBase->SetLineColor(colors[0]);
//   hBase->SetMarkerColor(colors[0]);
//   hBase->SetMarkerStyle(markers[0]);
//   hBase->SetFillColor(colors[0]);
//   hBase->SetMinimum(0.);
//   hBase->SetMaximum(1.);
//
//   canSuper->cd();
//   TH1* frame_super = canSuper->DrawFrame(0.,0.,8.,1.);
//   hBase->Draw("same");
//
//   canRatio->cd();
//   TH1* frame_ratio = canRatio->DrawFrame(0.,0.,8.,2.);
//   canBarlow->cd();
//   TH1* frame_barlow = canBarlow->DrawFrame(0.,0.,8.,10.);
//
//
//   TH1D* hTemp = 0x0;
//   TH1D* hRatio = 0x0;
//   // TH1D* hRatio2 = 0x0;
//   TH1D* hBarlow = 0x0;
//
//   for(Short_t iFile(1); iFile < iNumFiles; iFile++)
//   {
//     hTemp = (TH1D*) list->At(iFile);
//     if(!hTemp) { printf("ProcessList::Histo (%d) not found\n",iFile); return; }
//
//     hTemp->SetLineColor(colors[iFile]);
//     hTemp->SetMarkerColor(colors[iFile]);
//     hTemp->SetMarkerStyle(markers[iFile]);
//     hTemp->SetFillColor(colors[iFile]);
//
//     legendFull->AddEntry(hTemp,sTaskTag[iFile],"pl");
//     legend->AddEntry(hTemp,sTaskTag[iFile],"pl");
//
//     canSuper->cd();
//     hTemp->Draw("same");
//
//     hRatio = MakeRatio(hBase,hTemp,1);
//     // hRatio->SetLineColor(colors[iFile]);
//     // hRatio->SetMarkerColor(colors[iFile]);
//     // hRatio->SetMarkerStyle(markers[iFile]);
//     // hRatio->SetFillColor(colors[iFile]);
//
//     hBarlow = DoBarlowTest(hBase,hTemp,1);
//
//     canRatio->cd();
//     hRatio->Draw("hist p e1 same");
//
//     canBarlow->cd();
//     hBarlow->Draw("hist p same");
//   }
//
//   canSuper->cd();
//   legendFull->Draw();
//
//   // fitting the ratios
//   canRatio->cd();
//   legend->Draw();
//   TF1* fitConst = new TF1("fitConst","[0]",hRatio->GetXaxis()->GetXmin(),hRatio->GetXaxis()->GetXmax());
//   fitConst->SetLineColor(kRed);
//   hRatio->Fit("fitConst","R");
//   printf("Fit: %g ± %g (Chi2 %g)\n",fitConst->GetParameter(0),fitConst->GetParError(0),fitConst->GetChisquare());
//
//   TLine* lineUnity = new TLine();
//   lineUnity->SetLineColor(kBlack);
//   lineUnity->SetLineStyle(7);
//
//   fitConst->Draw("same");
//   lineUnity->DrawLine(hRatio->GetXaxis()->GetXmin(),1.,hRatio->GetXaxis()->GetXmax(),1.);
//
//   canBarlow->cd();
//   lineUnity->DrawLine(hBarlow->GetXaxis()->GetXmin(),1.,hBarlow->GetXaxis()->GetXmax(),1.);
//   legend->Draw();
//
//   return;
// }
//_____________________________________________________________________________
TH1D* MakeRatio(TH1D* hBase, TH1D* hSyst, Bool_t bCorelated)
{
  if(!hBase) { printf("MakeRatio::Baseline histo not found\n"); return 0x0; }
  if(!hSyst) { printf("MakeRatio::Systematics histo not found\n"); return 0x0; }

  TH1D* hRatio = (TH1D*) hSyst->Clone("hRatio");
  hRatio->Reset();

  // checking the binning
  Short_t iNumBins = hRatio->GetNbinsX();
  Short_t iNumBinsBase = hBase->GetNbinsX();
  Short_t iNumBinsSyst = hSyst->GetNbinsX();

  if( iNumBins != iNumBinsBase || iNumBins != iNumBinsSyst ) { printf("MakeRatio::Different number of bins\n"); return 0x0; }

  // calculating ratio
  Double_t dContBase = 0, dErrBase = 0;
  Double_t dContSyst = 0, dErrSyst = 0;
  Double_t dContRatio = 0, dErrRatio = 0;
  for(Short_t iBin(1); iBin < iNumBins+1; iBin++)
  {
    dContBase = hBase->GetBinContent(iBin);
    dErrBase = hBase->GetBinError(iBin);
    dContSyst = hSyst->GetBinContent(iBin);
    dErrSyst = hSyst->GetBinError(iBin);

    dContRatio =  dContSyst / dContBase;
    dErrRatio = TMath::Power(dErrSyst/dContBase, 2) + TMath::Power( dErrBase*dContSyst/(dContBase*dContBase), 2);
    // printf("Err (before) : %g | ", TMath::Sqrt(dErrRatio));

    if(bCorelated) dErrRatio -= (2*dContSyst*dErrBase*dErrSyst/TMath::Power(dContBase,3));
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    hRatio->SetBinContent(iBin,dContRatio);
    hRatio->SetBinError(iBin,TMath::Sqrt(dErrRatio));
  }

  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(1.5);

  return hRatio;
}
//_____________________________________________________________________________
TH1D* DoBarlowTest(TH1D* hBase, TH1D* hSyst, Bool_t bCorelated)
{
  if(!hBase) { printf("DoBarlowTest::Baseline histo not found\n"); return 0x0; }
  if(!hSyst) { printf("DoBarlowTest::Systematics histo not found\n"); return 0x0; }

  TH1D* hBarlow = (TH1D*) hSyst->Clone("hBarlow");
  hBarlow->Reset();

  // checking the binning
  Short_t iNumBins = hBarlow->GetNbinsX();
  Short_t iNumBinsBase = hBase->GetNbinsX();
  Short_t iNumBinsSyst = hSyst->GetNbinsX();

  if( iNumBins != iNumBinsBase || iNumBins != iNumBinsSyst ) { printf("DoBarlowTest::Different number of bins\n"); return 0x0; }

  // calculating ratio
  Double_t dContBase = 0, dErrBase = 0;
  Double_t dContSyst = 0, dErrSyst = 0;
  Double_t dBarlow = 0, dDenom = 0;
  for(Short_t iBin(1); iBin < iNumBins+1; iBin++)
  {
    dContBase = hBase->GetBinContent(iBin);
    dErrBase = hBase->GetBinError(iBin);
    dContSyst = hSyst->GetBinContent(iBin);
    dErrSyst = hSyst->GetBinError(iBin);

    if(bCorelated) { dDenom = TMath::Sqrt(TMath::Abs(TMath::Power(dErrBase,2) - TMath::Power(dErrSyst,2))); }
    else { dDenom = TMath::Sqrt(TMath::Abs(TMath::Power(dErrBase,2) + TMath::Power(dErrSyst,2))); }

    dBarlow = TMath::Abs(dContBase - dContSyst) / dDenom;

    hBarlow->SetBinContent(iBin,dBarlow);
    hBarlow->SetBinError(iBin,0);
  }

  // hBarlow->SetMinimum(0.);
  // hBarlow->SetMaximum(10.);

  return hBarlow;
}
