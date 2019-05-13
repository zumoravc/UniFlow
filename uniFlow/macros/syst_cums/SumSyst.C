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

Double_t dMax = 0.3;

void SetCustomPalette();
TH1D* MakeHist(TF1* fit, TH1D* histTemp);
void SumSystAll();

Bool_t ProcessSingle(
    const char* histNameC,
    const char* histNameLabel,
    const char* path,
    std::vector<TString>& vec
);

// ==================================================================================================================
// ==================================================================================================================

void SumSyst(TString species = "Lambda")
{
    Int_t iNumCent = 6;
    TString path = Form("/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/PbPb/syst_6815/%s/",species.Data());

    std::vector<TString> vecHist = {};
    vecHist.push_back("hFlow4_harm2_gap-10");

    std::vector<TString> vecHistLabel = {};
    vecHistLabel.push_back(Form("%s: v_{2}{4}",species.Data()));

    std::vector<TString> vecSyst = {};

    if(species.EqualTo("Charged")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("FB768");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");

    } else if(species.EqualTo("Pion")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("FB768");
        vecSyst.push_back("PID3sigma");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");

    } else if(species.EqualTo("Kaon")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("FB768");
        vecSyst.push_back("PID3sigma");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");

    } else if(species.EqualTo("Proton")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("FB768");
        vecSyst.push_back("PID3sigma");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");

    } else if(species.EqualTo("K0s")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");
        vecSyst.push_back("V0sCPA099");
        vecSyst.push_back("V0sCrossFind1");
        vecSyst.push_back("V0sDaugDCA3");
        vecSyst.push_back("V0sDaugPt02");
        vecSyst.push_back("V0sDecRad1");
        vecSyst.push_back("V0sDecRad10");
        vecSyst.push_back("V0sFinderOn");
        vecSyst.push_back("V0sPVDCA3");

    } else if(species.EqualTo("Lambda")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");
        vecSyst.push_back("V0sCPA099");
        vecSyst.push_back("V0sCrossFind1");
        vecSyst.push_back("V0sDaugDCA3");
        vecSyst.push_back("V0sDaugPt02");
        vecSyst.push_back("V0sDecRad1");
        vecSyst.push_back("V0sDecRad10");
        vecSyst.push_back("V0sFinderOn");
        vecSyst.push_back("V0sPVDCA3");

    } else if(species.EqualTo("Phi")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("FB768");
        vecSyst.push_back("PID3sigma");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");

        iNumCent = 5;
    } else {
        return;
    }

    for(Int_t iCent(0); iCent < iNumCent; ++iCent) {
        // Int_t iCent = 2;
        for(Int_t iHist(0); iHist < (Int_t) vecHist.size(); ++iHist) {
            // Int_t iHist = 0;

            TString histoName = vecHist.at(iHist);
            TString histoNameLabel = vecHistLabel.at(iHist);

            if(!ProcessSingle(Form("%s_%s_cent%d",species.Data(),histoName.Data(),iCent),histoNameLabel.Data(),path.Data(),vecSyst)) { return; }
        }
    }

    return;
}

// ==================================================================================================================
// ==================================================================================================================

Bool_t ProcessSingle(const char* histNameC,const char* histNameLabel, const char* path, std::vector<TString>& vec)
{
    TString histName = histNameC;

    Int_t iNumSyst = vec.size();
    if(iNumSyst == 0) { printf("E: Input vector is empty!\n"); return kFALSE; }

    std::vector<TF1*> vecObj;
    std::vector<TFile*> vecFiles;

    // template (ranges) histo
    TH1D* hist = nullptr;

    for(Int_t iSyst(0); iSyst < iNumSyst; ++iSyst) {

        TString sSyst = vec.at(iSyst);
        TString sObjectName = "fitDiff";

        TFile* fileIn = TFile::Open(Form("%s/%s/syst_root/%s.root",path,sSyst.Data(),histName.Data()),"READ");
        if(!fileIn) { printf("E: Input file '%s' not found!\n", histName.Data()); return kFALSE; }

        TList* list = (TList*) fileIn->Get("list");
        if(!list) { printf("E: Input list not found!\n"); fileIn->ls(); return kFALSE; }

        TF1* fit = (TF1*) list->FindObject(sObjectName.Data());
        if(!fit) { printf("E: Input object '%s' not found!\n",sObjectName.Data()); list->ls(); return kFALSE; }

        vecObj.push_back(fit);
        vecFiles.push_back(fileIn);

        if(!hist) {
            hist = (TH1D*) list->FindObject("Diff");
            if(!hist) {
                printf("E: Template input histo not found!\n");
                list->ls();
                return kFALSE;
            }
        }
    }

    Int_t iNumObj = vecObj.size();
    printf("Found %d objects!\n", iNumObj);
    if(iNumObj != iNumSyst) { printf("E: Unexpected number of objects!\n"); return kFALSE; }

    // plotting stuff

    hist->Reset();
    hist->SetTitle(histNameLabel);
    hist->GetYaxis()->SetTitle("Syst. uncrt. (%)");
    hist->SetMaximum(dMax);
    hist->SetLabelFont(43,"XY");
    hist->SetLabelSize(27,"XY");
    hist->SetTitleFont(43,"XY");
    hist->SetTitleSize(27,"XY");
    hist->SetTitleOffset(1.2,"Y");
    hist->SetTitleOffset(1.2,"X");

    TLegend* leg = new TLegend(0.74,0.16,0.88,0.88);
    leg->SetTextFont(43);
    leg->SetTextSize(18);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0,0);
    // leg->AddEntry(histBase,"Baseline","pel");
    // leg->AddEntry(histSyst,syst,"pel");

    TLine* lineUnity = new TLine();
    lineUnity->SetLineStyle(2);
    lineUnity->SetLineWidth(3);
    lineUnity->SetLineColor(kGray+1);

    TF1* fitTotal = (TF1*) vecObj.at(0)->Clone(Form("%s",histName.Data()));
    fitTotal->SetLineColor(kRed);
    fitTotal->SetLineStyle(kSolid);
    fitTotal->SetLineWidth(4);
    leg->AddEntry(fitTotal,"Total","l");

    SetCustomPalette();
    Int_t nnCol = gStyle->GetNumberOfColors();
    Int_t nPnt  = iNumSyst;

    TCanvas* can = new TCanvas("can","can",1000,700);
    can->cd(1);
    gPad->SetRightMargin(0.3);
    hist->Draw("hist");

    Double_t dTotalSq = 0.0;

    for(Int_t i(0); i < iNumObj; ++i) {

        Int_t idx = i * Float_t(nnCol-1) / (nPnt-1);

        TF1* fit = vecObj.at(i);
        fit->SetLineColor(gStyle->GetColorPalette(idx));
        fit->SetLineStyle(kSolid);
        fit->SetLineWidth(3);
        fit->Draw("same");

        Double_t dSyst = fit->GetParameter(0);
        dTotalSq += dSyst*dSyst;

        leg->AddEntry(fit,vec.at(i).Data(),"l");
    }

    fitTotal->SetParameter(0,TMath::Sqrt(dTotalSq));

    fitTotal->Draw("same");
    lineUnity->DrawLine(hist->GetXaxis()->GetXmin(), 0.0, hist->GetXaxis()->GetXmax(), 0.0);
    leg->Draw();

    TH1D* histTotal = MakeHist(fitTotal, hist);
    if(!histTotal) {printf("E: histTotal not exits!\n"); return kFALSE; }

    gSystem->mkdir(Form("%s/final_fit/syst_pdf/",path),kTRUE);
    // gSystem->mkdir(Form("%s/final_fit/syst_eps/",path),kTRUE);
    can->SaveAs(Form("%s/final_fit/syst_pdf/%s.pdf",path,histName.Data()),"pdf");
    // can->SaveAs(Form("%s/final_fit/syst_eps/%s.eps",path,histName.Data()),"eps");

    TFile* fileOut = TFile::Open(Form("%s/final_fit/syst_total.root",path),"UPDATE");
    if(!fileOut) { printf("E: Output file not opened!\n"); return kFALSE; }
    fileOut->cd();
    // fitTotal->Write(Form("%s",histName.Data()));
    histTotal->Write(Form("%s",histName.Data()));
    fileOut->Close();

    for(Int_t iFile(0); iFile < vecFiles.size(); ++iFile) {
        TFile* f = vecFiles.at(iFile);
        f->Close();
    }


    return kTRUE;
}
// ==================================================================================================================
TH1D* MakeHist(TF1* fit, TH1D* histTemp)
{
    if(!fit) { printf("ERROR-MakeHist: no fit!\n"); return nullptr; }
    if(!histTemp) { printf("ERROR-MakeHist: no histTemp!\n"); return nullptr; }

    Double_t dC = fit->GetParameter(0);

    TH1D* histTotal = (TH1D*) histTemp->Clone("histTotal");
    for(Int_t iBin(1); iBin < histTemp->GetNbinsX()+1; ++iBin) {
        histTotal->SetBinContent(iBin, dC);
    }

    return histTotal;
}
// ==================================================================================================================
void SumSystAll()
{
    SumSyst("Charged");
    SumSyst("Pion");
    SumSyst("Kaon");
    SumSyst("Proton");
    SumSyst("K0s");
    SumSyst("Lambda");
    SumSyst("Phi");
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
