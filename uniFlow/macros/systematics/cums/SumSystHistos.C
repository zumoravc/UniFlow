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

#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/utils/Utils.cxx>

std::vector<Double_t> vecPID = {};
std::vector<Double_t> vecDec = {};

Double_t dMax = 0.2;
Bool_t AddHistos(std::vector<TH1D*>& vec, TH1D* total);
TH1D* MakeValues(Double_t dval, TH1D* temp);
void  FitToHist(TH1* hist, TF1* fit);
void SetCustomPalette();

Bool_t ProcessSingle(
    const char* histNameC,
    const char* path,
    std::vector<TString>& vec,
    Int_t iCent,
    TString sSystHistoName = "Diff",
    TString sSystFitName = ""
);

void ProcessAll();

// ==================================================================================================================
// ==================================================================================================================

void SumSystHistos(
    TString species = "Lambda",
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst/",
    Int_t iNumCent = 6,
    TString sSystHistoName = "Diff",
    TString sSystFitName = "fitDiff"
)
{
    TString path = Form("%s/%s/",sPath.Data(),species.Data());

    std::vector<TString> vecHist = {"hFlow4_harm2_gap-10"};
    // vecHist.push_back("<<3>>(4,-2,-2)_2sub(0)");
    // vecHist.push_back("<<3>>(5,-3,-2)_2sub(0)");
    // vecHist.push_back("<<3>>(6,-3,-3)_2sub(0)");

    std::vector<TString> vecSyst;
    if(species.EqualTo("Charged")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("FB768");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");

    } else if(species.EqualTo("Pion")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("FB768");
        vecSyst.push_back("PID2sigma");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");
        // vecPID = {0.01,0.02,0.02,0.04,0.03,0.04};

    } else if(species.EqualTo("Kaon")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("FB768");
        vecSyst.push_back("PID2sigma");
        vecSyst.push_back("Bayes90");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");
        // vecPID = {0.04,0.02,0.01,0.01,0.01,0.01};

    } else if(species.EqualTo("Proton")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("FB768");
        vecSyst.push_back("PID2sigma_anti");
        vecSyst.push_back("Bayes90");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");
        // vecPID = {0.02,0.03,0.03,0.03,0.03,0.03};

    } else if(species.EqualTo("K0s")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");
        // vecSyst.push_back("V0sCPA099");
        vecSyst.push_back("V0sCrossFind1");
        vecSyst.push_back("V0sDaugDCA3");
        vecSyst.push_back("V0sDaugPt02");
        vecSyst.push_back("V0sDecRad1");
        vecSyst.push_back("V0sDecRad10");
        vecSyst.push_back("V0sFinderOn");
        vecSyst.push_back("V0sPVDCA3");
        // vecDec = {0.02,0.024,0.02,0.012,0.021,0.036};

    } else if(species.EqualTo("Lambda")) {
        vecSyst.push_back("CL1");
        vecSyst.push_back("PVz8");
        vecSyst.push_back("TPCcls90");
        // vecSyst.push_back("V0sCPA099");
        vecSyst.push_back("V0sCrossFind1");
        vecSyst.push_back("V0sDaugDCA3");
        vecSyst.push_back("V0sDaugPt02");
        vecSyst.push_back("V0sDecRad1");
        vecSyst.push_back("V0sDecRad10");
        vecSyst.push_back("V0sFinderOn");
        vecSyst.push_back("V0sPVDCA3");

    } else if(species.EqualTo("Phi")) {
        vecSyst.push_back("CL1");
        // vecSyst.push_back("FB768");
        // vecSyst.push_back("PID2sigma");
        // vecSyst.push_back("Bayes90");
        // vecSyst.push_back("PVz8");
        // vecSyst.push_back("TPCcls90");
        // vecPID = {0.04,0.02,0.01,0.01,0.01,0.01};

        // iNumCent = 6;
    } else {
        return;
    }

    printf("iNumCent %d\n",iNumCent);
    for(Int_t iCent(0); iCent < iNumCent; ++iCent) {
        printf("iCent %d\n",iCent);
        for(Int_t iHist(0); iHist < (Int_t) vecHist.size(); ++iHist) {
            // Int_t iHist = 0;

            TString histoName = vecHist.at(iHist);

            if(!ProcessSingle(Form("%s_%s_cent%d",species.Data(),histoName.Data(),iCent),path.Data(),vecSyst,iCent,sSystHistoName,sSystFitName)) { return; }
        }
    }

    return;

}

// ==================================================================================================================
// ==================================================================================================================

Bool_t ProcessSingle(
    const char* histNameC,
    const char* path,
    std::vector<TString>& vec,
    Int_t iCent,
    TString sSystHistoName,
    TString sSystFitName
)
{
    Bool_t bUseSystFit = sSystFitName.Length();

    TString histName = histNameC;

    Int_t iNumSyst = vec.size();
    if(iNumSyst == 0) { Utils::Debug("Input vector is empty!"); return kFALSE; }

    TString sOutDir = "final_combined";

    std::vector<TH1D*> vecObj;
    std::vector<TFile*> vecFiles;

    // template (ranges) histo
    TH1D* hist = nullptr;
    TF1* fit = nullptr;

    for(Int_t iSyst(0); iSyst < iNumSyst; ++iSyst) {

        TString sSyst = vec.at(iSyst);
        Utils::Debug(Form("Adding Syst Source %s",sSyst.Data()));

        TFile* fileIn = TFile::Open(Form("%s/%s/syst_root/%s.root",path,sSyst.Data(),histName.Data()),"READ");
        if(!fileIn) { Utils::Error(Form("Input file '%s' not found!", histName.Data())); return kFALSE; }
        vecFiles.push_back(fileIn);

        TList* list = (TList*) fileIn->Get("list");
        if(!list) { Utils::Error("Input list 'list' not found!"); fileIn->ls(); return kFALSE; }

        hist = (TH1D*) list->FindObject(sSystHistoName.Data());

        if(sSyst.EqualTo("PID3sigma") && vecPID.size() > 0) {
            hist = MakeValues(vecPID.at(iCent),hist);
            bUseSystFit = kFALSE;
        }

        if(sSyst.EqualTo("V0sDecRad10") && vecDec.size() > 0) {
            hist = MakeValues(vecDec.at(iCent),hist);
            bUseSystFit = kFALSE;
        }

        if(!hist) {
            Utils::Error("Syst histo not found!");
            list->ls();
            return kFALSE;
        }

        vecObj.push_back(hist);

        if(bUseSystFit) {
            fit = (TF1*) list->FindObject(sSystFitName.Data());
            if(!fit) { Utils::Error(Form("Syst fit '%s' not found!",sSystFitName.Data())); list->ls(); return kFALSE; }
            FitToHist(hist,fit);
        }
    }

    Int_t iNumObj = vecObj.size();
    Utils::Info(Form("Found %d objects!", iNumObj));
    if(iNumObj != iNumSyst) { Utils::Error("Unexpected number of objects!"); return kFALSE; }

    // plotting stuff

    TLegend* leg = new TLegend(0.74,0.16,0.88,0.88);
    leg->SetTextFont(43);
    leg->SetTextSize(18);
    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0,0);

    TLine* lineUnity = new TLine();
    lineUnity->SetLineStyle(2);
    lineUnity->SetLineWidth(2);
    lineUnity->SetLineColor(kGray+1);

    TH1D* fitTotal = (TH1D*) vecObj.at(0)->Clone(Form("%s_total",histName.Data()));
    fitTotal->Reset();
    if(!AddHistos(vecObj,fitTotal)) return kFALSE;
    fitTotal->SetLineColor(kRed);
    fitTotal->SetLineStyle(kDashed);
    fitTotal->SetLineWidth(3);
    leg->AddEntry(fitTotal,"Total","l");

    SetCustomPalette();
    Int_t nnCol = gStyle->GetNumberOfColors();
    Int_t nPnt  = iNumSyst;

    TCanvas* can = new TCanvas("can","can",1000,700);
    can->cd(1);
    gPad->SetRightMargin(0.3);
    TH1* frame = (TH1*) gPad->DrawFrame(hist->GetXaxis()->GetXmin(),0.0,hist->GetXaxis()->GetXmax(),dMax);
    // hist->Reset();
    frame->SetTitle(histName.Data());
    // frame->SetMaximum(dMax);
    frame->SetLabelFont(43,"XY");
    frame->SetLabelSize(27,"XY");
    frame->SetTitleFont(43,"XY");
    frame->SetTitleSize(27,"XY");
    frame->SetTitleOffset(1.2,"Y");
    frame->SetTitleOffset(1.2,"X");

    for(Int_t i(0); i < iNumObj; ++i) {

        Int_t idx = i * Float_t(nnCol-1) / (nPnt-1);

        TH1D* hist = vecObj.at(i);
        hist->SetLineColor(gStyle->GetColorPalette(idx));
        hist->SetLineWidth(3);
        hist->SetLineStyle(kSolid);
        hist->Draw("hist same");

        // Double_t dSyst = fit->GetParameter(0);

        leg->AddEntry(hist,vec.at(i).Data(),"l");
    }


    fitTotal->Draw("same");
    // lineUnity->DrawLine(hist->GetXaxis()->GetXmin(), 0.0, hist->GetXaxis()->GetXmax(), 0.0);
    leg->Draw();

    gSystem->mkdir(Form("%s/%s/syst_pdf/",path,sOutDir.Data()),kTRUE);
    // gSystem->mkdir(Form("%s/final/syst_eps/",path),kTRUE);
    can->SaveAs(Form("%s/%s/syst_pdf/%s.pdf",path,sOutDir.Data(), histName.Data()),"pdf");
    // can->SaveAs(Form("%s/final/syst_eps/%s.eps",path,histName.Data()),"eps");

    TFile* fileOut = TFile::Open(Form("%s/%s/syst_total.root",path,sOutDir.Data()),"UPDATE");
    if(!fileOut) { printf("E: Output file not opened!\n"); return kFALSE; }
    fileOut->cd();
    fitTotal->Write(Form("%s",histName.Data()));
    fileOut->Close();

    for(Int_t iFile(0); iFile < vecFiles.size(); ++iFile) {
        TFile* f = vecFiles.at(iFile);
        f->Close();
    }


    return kTRUE;
}
// ==================================================================================================================
Bool_t AddHistos(std::vector<TH1D*>& vec, TH1D* total)
{
    if(!total) { Utils::Error("hist total not found!","AddHistos"); return kFALSE; }

    Int_t iNumBins = total->GetNbinsX();
    Int_t iNumHist = vec.size();

    Utils::Debug(Form("Adding %d",iNumHist),1,"AddHistos");

    for(Int_t i(1); i < iNumBins+1; ++i) {
        Double_t dValue = 0.0;

        for(Int_t h(0); h < iNumHist; ++h) {
            TH1D* hist = (TH1D*) vec.at(h);
            if(!hist) { Utils::Error(Form("Hist %d not found!",h)); return kFALSE; }
            Double_t d = hist->GetBinContent(i);
            Utils::Debug(Form("Syst source %d | bin %d | Value %g",h,i,d),1,"AddHistos");

            dValue += d*d;
        }

        dValue = TMath::Sqrt(dValue);
        total->SetBinContent(i, dValue);
        Utils::Debug(Form("SystTotal | bin %d | Value %g",i,dValue),1,"AddHistos");
    }

    return kTRUE;
}
// ==================================================================================================================
TH1D* MakeValues(Double_t val, TH1D* temp)
{
    if(!temp) { printf("MakeValues : histo not found!\n"); return nullptr; }
    TH1D* hist = (TH1D*) temp->Clone("histTemp");
    hist->Reset();

    printf("%g\n",val);

    for(Int_t iBin(1); iBin < temp->GetNbinsX()+1; ++iBin) {
        hist->SetBinContent(iBin,val);
    }

    return hist;
}
// ==================================================================================================================
void  FitToHist(TH1* hist, TF1* fit)
{
    if(!hist) { Utils::Error("Input histo not found!","FitToHist"); return; }
    if(!fit) { Utils::Error("Input fit not found!"); return; }

    Double_t dFitXmin = fit->GetXmin();
    Double_t dHistXmin = hist->GetXaxis()->GetXmin();

    Double_t dFitValue = fit->GetParameter(0);

    Int_t iBins = hist->GetNbinsX();

    Bool_t bIsAboveCut = kFALSE;
    for(Int_t iB(1); iB < iBins+1; ++iB) {
        if(!bIsAboveCut) {
            Double_t dEdgeLow = hist->GetBinLowEdge(iB);
            Double_t dEdgeHigh = hist->GetBinLowEdge(iB+1);

            if(dEdgeHigh > dFitXmin) {
                bIsAboveCut = kTRUE;
                Utils::Debug(Form("First bin above cut %g (bin %g-%g)",dFitXmin,dEdgeLow,dEdgeHigh),1,"FitToHist");
            }
        }

        if(bIsAboveCut) {
            hist->SetBinContent(iB, dFitValue);
        }
    }

    return;
}
// ==================================================================================================================
void ProcessAll()
{
    SumSystHistos("Charged");
    SumSystHistos("Pion");
    SumSystHistos("Kaon");
    SumSystHistos("Proton");
    SumSystHistos("K0s");
    SumSystHistos("Lambda");
    SumSystHistos("Phi");
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
