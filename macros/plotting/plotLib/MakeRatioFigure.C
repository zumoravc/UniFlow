#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TGaxis.h"

#include "TH1.h"

// Example displaying two histograms and their ratio.
// Based on ROOT tutorial macro by Olivier Couet

void PrepareCanvas(TPad* padEnvelope, TPad* padMain, TPad* padRatio);

void MakeRatioFigure()
{

  // Define two gaussian histograms. Note the X and Y title are defined
  // at booking time using the convention "Hist_title ; X_title ; Y_title"
  TH1F *h1 = new TH1F("h1", "Two gaussian plots and their ratio;x title; h1 and h2 gaussian histograms", 100, -5, 5);
  TH1F *h2 = new TH1F("h2", "h2", 100, -5, 5);
  h1->FillRandom("gaus");
  h2->FillRandom("gaus");


  // NEW IMPLEMENTATION ==================================================
  TCanvas* canNew = new TCanvas("canNew","canNew",800,800);
  canNew->Divide(1,1);
  // canNew->cd(1);

  TPad* padMain = 0x0;
  TPad* padRatio = 0x0;

  PrepareCanvas((TPad*)canNew->cd(1),padMain,padRatio);

  canNew->cd();
  // padMain->cd(1);
  // gPad->SetFillColor(kGreen);
  h1->Draw("same");
  // padMain->SetFillColor(kAzure);

}

void PrepareCanvas(TPad* padEnvelope, TPad* padMain, TPad* padRatio)
{
  padEnvelope->cd();


  padMain = new TPad("padMain","padMain", 0, 0.3, 1, 1.0);
  padMain->SetBottomMargin(0.015);
  padMain->SetFillColor(kRed);
  padMain->Draw();
  padMain->cd();
  TH1* frame_main = (TH1*) gPad->DrawFrame(1.0,0.0,10.0,220.0);
  frame_main->SetTitle("Title; x; y");
  TAxis* axis_x_main = frame_main->GetXaxis();
  TAxis* axis_y_main = frame_main->GetYaxis();
  axis_y_main->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis_y_main->SetLabelSize(15);
  axis_y_main->SetNdivisions(510);
  axis_x_main->SetLabelSize(0.);

  axis_y_main->SetTitleSize(20);
  axis_y_main->SetTitleFont(43);
  axis_y_main->SetTitleOffset(1.55);

  padEnvelope->cd();
  padRatio = new TPad("padRatio","padRatio", 0, 0.0, 1, 0.3);
  padRatio->SetTopMargin(0.025);
  padRatio->SetBottomMargin(0.25);
  padRatio->SetFillColor(kBlue-9);
  padRatio->Draw();
  padRatio->cd();
  TH1* frame_ratio = (TH1*) gPad->DrawFrame(1.0,0.0,10.0,2.0);
  frame_ratio->SetTitle("; x; y");
  TAxis* axis_x_ratio = frame_ratio->GetXaxis();
  TAxis* axis_y_ratio = frame_ratio->GetYaxis();
  // axis_y_main->SetLabelSize(0.0);
  axis_y_ratio->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis_y_ratio->SetLabelSize(15);
  axis_y_ratio->SetNdivisions(5,5,0,kTRUE);
  axis_x_ratio->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis_x_ratio->SetLabelSize(15);
  // axis_x_ratio->SetNdivisions(510);





  // Y axis ratio plot settings
  axis_y_ratio->SetTitle("ratio h1/h2 ");
  axis_y_ratio->SetNdivisions(505);
  axis_y_ratio->SetTitleSize(20);
  axis_y_ratio->SetTitleFont(43);
  axis_y_ratio->SetTitleOffset(1.55);
  axis_y_ratio->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis_y_ratio->SetLabelSize(15);

  // X axis ratio plot settings
  axis_x_ratio->SetTitleSize(20);
  axis_x_ratio->SetTitleFont(43);
  axis_x_ratio->SetTitleOffset(4.);
  axis_x_ratio->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axis_x_ratio->SetLabelSize(15);
}
