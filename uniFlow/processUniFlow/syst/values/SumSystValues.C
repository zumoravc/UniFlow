#include <vector>
#include "TSystem.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TText.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TF1.h"

void SetCustomPalette();

// TString corName = "v422";
// TString corName = "v523";
// TString corName = "v633";


// TString species = "K0s";
// TString species = "Lambda";
// TString species = "Phi";


void SumSystValues()
{

    TString species = "Lambda"; Double_t xMin = 0.8; Double_t xMax = 6.0;

    TString sSyst[] = {
        // "FB768",
        // "PID3sigma",
        "PVz8",
        "TPCcls90",
        "V0sCPA099",
        "V0sCrossFind1",
        "V0sDaugDCA3",
        "V0sDaugPt02",
        "V0sDecRad10",
        "V0sFinderOn",
        "V0sPVDCA3"
    };

    // TString corName = "v422";
    // Double_t dVal[] = {  0,     2,      0,      0,      1,      1,      2,      2,      2 };  Int_t iCent = 0; TString sCent = "0-5%";
    // Double_t dVal[] = {  2,     2,      0,      1,      2,      1,      2.1,    3,      3 };  Int_t iCent = 1; TString sCent = "5-10%";
    // Double_t dVal[] = {  1.2,   2.2,      2,      3,      1,      1,      3.2,   3,      4 };  Int_t iCent = 2; TString sCent = "10-20%";
    // Double_t dVal[] = {  1.1,   2,      1,      2,      1,      1,      3,      3,      4 };  Int_t iCent = 3; TString sCent = "20-30%";
    // Double_t dVal[] = {  1,     1,      1,      2,      1,      1,      2,      2,      2 };  Int_t iCent = 4; TString sCent = "30-40%";
    // Double_t dVal[] = {  0.1,   2.1,      0,      0.5,    1,      0.5,    1,      2,      2 };  Int_t iCent = 5; TString sCent = "40-50%";
    // Double_t dVal[] = {  2,     2,      0,      1,      2,      0.5,    2,      3,      3 };  Int_t iCent = 6; TString sCent = "50-60%";

    TString corName = "v523";
    // Double_t dVal[] = {  0,  2,  0,  1,  0,  2,  2,  1,  2 };  Int_t iCent = 0; TString sCent = "0-5%";
    // Double_t dVal[] = {  1,  2,  3,  1,  3,  2,  2,  1,  2.2 };  Int_t iCent = 1; TString sCent = "5-10%";
    // Double_t dVal[] = {  2,  2,  1,  2,  4,  3,  1,  1,  2.4 };  Int_t iCent = 2; TString sCent = "10-20%";
    // Double_t dVal[] = {  2,  2,  1,  2,  3,  3,  2,  1,  3 };  Int_t iCent = 3; TString sCent = "20-30%";
    // Double_t dVal[] = {  3,  2,  2,  1.4,  2,  2,  1,  1,  2.1 };  Int_t iCent = 4; TString sCent = "30-40%";
    // Double_t dVal[] = {  1,  2,  0,  1.1,  4,  3,  2,  1,  2 };  Int_t iCent = 5; TString sCent = "40-50%";
    Double_t dVal[] = {  0,  2,  0,  1,  0,  2,  0,  1,  2 };  Int_t iCent = 6; TString sCent = "50-60%";
    //
    // TString corName = "v633";
    // Double_t dVal[] = {  1,  2,  1,  3,  0,  0,  2,  1,  2 };  Int_t iCent = 0; TString sCent = "0-5%";
    // Double_t dVal[] = {  1,  2,  1,  3,  3,  3,  2,  1,  2 };  Int_t iCent = 1; TString sCent = "5-10%";
    // Double_t dVal[] = {  1,  2,  1,  3,  4,  3,  2,  1,  2 };  Int_t iCent = 2; TString sCent = "10-20%";
    // Double_t dVal[] = {  2,  2,  1,  3,  4,  2,  2,  1,  3 };  Int_t iCent = 3; TString sCent = "20-30%";
    // Double_t dVal[] = {  3,  2,  1,  3,  4,  3,  2,  1,  2 };  Int_t iCent = 4; TString sCent = "30-40%";
    // Double_t dVal[] = {  1,  2,  1,  3,  0,  0,  2,  1,  2 };  Int_t iCent = 5; TString sCent = "40-50%";
    // Double_t dVal[] = {  1,  2,  1,  3,  0,  0,  2,  1,  2 };  Int_t iCent = 6; TString sCent = "50-60%";


















    TString path = Form("/mnt/CodesALICE/Flow/uniFlow/results/nlf/systematics/%s/",species.Data());
    TString histName = Form("%s_%s",species.Data(),corName.Data());
    TString histTitle = Form("%s: %s (%s)",species.Data(),corName.Data(),sCent.Data());

    Int_t iNumSyst = sizeof(sSyst)/sizeof(sSyst[0]);

    TLegend* leg = new TLegend(0.74,0.16,0.88,0.88);
    leg->SetTextFont(43);
    leg->SetTextSize(18);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0,0);

    SetCustomPalette();
    Int_t nnCol = gStyle->GetNumberOfColors();
    Int_t nPnt  = iNumSyst;

    TCanvas* can = new TCanvas("can","can",1000,700);
    can->cd(1);
    gPad->SetRightMargin(0.3);
    TH1* hist = (TH1*) gPad->DrawFrame(xMin,0.0,xMax,0.3);
    // hist->Reset();
    hist->SetTitle(histTitle.Data());
    hist->SetLabelFont(43,"XY");
    hist->SetLabelSize(27,"XY");
    hist->SetTitleFont(43,"XY");
    hist->SetTitleSize(27,"XY");
    hist->SetTitleOffset(1.2,"Y");
    hist->SetTitleOffset(1.2,"X");

    TList* list = new TList();

    Double_t dTotalSq = 0.0;

    for(Int_t i(0); i < iNumSyst; ++i) {

        Int_t idx = i * Float_t(nnCol-1) / (nPnt-1);

        TString sSource = sSyst[i];
        Double_t dValue = 1.0 * dVal[i] / 100.0;

        TF1* fit = new TF1(sSource.Data(),"[0]",xMin,xMax);
        fit->FixParameter(0, dValue);
        fit->SetLineColor(gStyle->GetColorPalette(idx));
        fit->SetLineStyle(kSolid);
        fit->Draw("same");

        dTotalSq += dValue*dValue;

        leg->AddEntry(fit,sSource.Data(),"l");

        list->Add(fit);
    }

    TF1* fitTotal = new TF1("total","[0]",xMin,xMax);
    fitTotal->SetLineColor(kRed);
    fitTotal->SetLineStyle(kSolid);
    fitTotal->SetLineWidth(3);
    fitTotal->SetParameter(0,TMath::Sqrt(dTotalSq));
    leg->AddEntry(fitTotal,"Total","l");
    list->Add(fitTotal);
    // fitTotal->Write(Form("%s",histName.Data()));

    TLine* lineUnity = new TLine();
    lineUnity->SetLineStyle(2);
    lineUnity->SetLineWidth(2);
    lineUnity->SetLineColor(kGray+1);

    fitTotal->Draw("same");
    lineUnity->DrawLine(hist->GetXaxis()->GetXmin(), 0.0, hist->GetXaxis()->GetXmax(), 0.0);
    leg->Draw();

    gSystem->mkdir(Form("%s/final/values/syst_pdf/",path.Data()),kTRUE);
    gSystem->mkdir(Form("%s/final/values/syst_eps/",path.Data()),kTRUE);
    can->SaveAs(Form("%s/final/values/syst_pdf/%s_cent%d.pdf",path.Data(),histName.Data(),iCent),"pdf");
    can->SaveAs(Form("%s/final/values/syst_eps/%s_cent%d.eps",path.Data(),histName.Data(),iCent),"eps");


    TFile* fileOut = TFile::Open(Form("%s/final/values/%s.root",path.Data(),histName.Data()),"UPDATE");
    if(!fileOut) { printf("E: Output file not opened!\n"); return; }
    fileOut->cd();
    list->Write(Form("mult%d",iCent),TObject::kSingleKey+TObject::kOverwrite);

    return;
}
// ==================================================================================================================
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
