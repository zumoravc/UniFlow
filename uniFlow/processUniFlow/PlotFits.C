// Macro for plotting fits within fits.root file from ProcessUniFlow
// March 2018 (for NLF modes systematics run)

#include <vector>
#include "TCanvas.h"
#include "TSystem.h"
#include "TPad.h"
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"

TCanvas* gCanSingle = nullptr;
TCanvas* gCanOverMass = nullptr;
TCanvas* gCanOverCorr = nullptr;

Double_t markSize = 1.1;
Int_t lineSize = 2.0;

Color_t colorTot = kRed;
Color_t colorSig = kGreen+2;
Color_t colorBg = kBlue;

Style_t styleTot = kSolid;
Style_t styleSig = kDashed;
Style_t styleBg = kDashed;

TString sFileName = "fits.root";
TString gOutFormat = "pdf";

// Double_t dMax = 1e-4;
// Double_t dMin = -1.0*1e-4;

void PrepareCanvas(Int_t iNumPt);
void SetPad();
void SetHistOff(TH1* hist);
void SetFuncAtt(TF1* func, Color_t color, Style_t style, Int_t width = 2.0);
void SetHistAtt(TH1* hist, Color_t color, Style_t markStyle, Double_t markSize = 1.0);
Bool_t ProcessList(TFile* file, TString sPath, TString sListName, Int_t padIndex = -1);

void PlotFits(TString sPath = "../results/nlf/output/K0s/", TString corrName = "<<2>>(2,-2)", TString sSpecies = "K0s", Int_t iNumPt = 8, Int_t iNumCent = 5)
{

    PrepareCanvas(iNumPt);


    TFile* fInput = TFile::Open(Form("%s/%s",sPath.Data(),sFileName.Data()),"READ");
    if(!fInput) { printf("ERROR: File not open!\n"); return; }

    gSystem->mkdir(Form("%s/fits/",sPath.Data()),1);

    TString sCorrName = corrName;

    for(Int_t iCent(0); iCent < iNumCent; ++iCent) {

        for(Int_t iPt(0); iPt < iNumPt; ++iPt) {

            TString name = Form("%s_%s_cent%d_pt%d",sSpecies.Data(),sCorrName.Data(),iCent,iPt);
            if(!ProcessList(fInput, sPath, name, iPt+1)) { printf("ERROR: ProcessList '%s' failed!\n", name.Data()); return; }
        }

        gCanOverMass->cd(0);
        gCanOverMass->SaveAs(Form("%s/fits/mass/mass_%s_%s_cent%d.%s",sPath.Data(),sSpecies.Data(),sCorrName.Data(),iCent,gOutFormat.Data()),gOutFormat.Data());
        gCanOverMass->Clear("D");
        gCanOverCorr->cd(0);
        gCanOverCorr->SaveAs(Form("%s/fits/corr/corr_%s_%s_cent%d.%s",sPath.Data(),sSpecies.Data(),sCorrName.Data(),iCent,gOutFormat.Data()),gOutFormat.Data());
        gCanOverCorr->Clear("D");

    }

    // if(!ProcessList(fInput, listName)) { printf("ERROR: ProcessList '%s' failed!\n", listName.Data()); return; }
    return;
}


Bool_t ProcessList(TFile* file, TString sPath, TString sListName, Int_t padIndex)
{
    if(!file) { printf("ERROR: Input file not found!\n"); return kFALSE; }

    TList* list = (TList*) file->Get(sListName.Data());
    if(!list) { printf("ERROR: Input list '%s' not found!\n",sListName.Data()); file->ls(); return kFALSE; }

    TString sCorr = ((TNamed*) list->FindObject("corr"))->GetTitle();
    TString sSpec = ((TNamed*) list->FindObject("spec"))->GetTitle();
    TString sCent = ((TNamed*) list->FindObject("cent"))->GetTitle();
    TString sPt = ((TNamed*) list->FindObject("pt"))->GetTitle();

    TString sTitle = Form("%s: %s GeV (%s%%)",sSpec.Data(), sPt.Data(), sCent.Data());

    // ========== Processing histCorr =======================================

    TH1D* histCorr = (TH1D*) list->FindObject("histCorr");
    if(!histCorr) { printf("ERROR: Object '%s' not found!\n", "histCorr"); list->ls(); return kFALSE; }

    TF1* fitCorr = (TF1*) list->FindObject("fitCorr");
    if(!fitCorr) { printf("ERROR: Object '%s' not found!\n", "fitCorr"); list->ls(); return kFALSE; }

    TF1* fitCorrSig = (TF1*) list->FindObject("fitCorSig");
    if(!fitCorrSig) { printf("ERROR: Object '%s' not found!\n", "fitCorrSig"); list->ls(); return kFALSE; }

    TF1* fitCorrBg = (TF1*) list->FindObject("fitCorBg");
    if(!fitCorrBg) { printf("ERROR: Object '%s' not found!\n", "fitCorrBg"); list->ls(); return kFALSE; }

    Double_t dXlow = fitCorr->GetXmin();
    Double_t dXhigh = fitCorr->GetXmax();

    histCorr->SetTitle(sTitle.Data());
    histCorr->GetXaxis()->SetRangeUser(dXlow,dXhigh);
    histCorr->SetStats(0);
    // histCorr->SetMinimum(0.0);

    // histCorr->SetMinimum(-0.05*histCorr->GetBinContent(histCorr->GetMaximumBin()));
    // histCorr->SetMaximum(1.2*histCorr->GetBinContent(histCorr->GetMaximumBin()));

    histCorr->SetMinimum();
    histCorr->SetMaximum();
    if(histCorr->GetMinimum() > 0.0) { histCorr->SetMinimum(0.0); }


    // histCorr->SetMaximum(dMax);
    // histCorr->SetMaximum(0.0);
    SetHistAtt(histCorr, kBlack, kDot, markSize);

    SetFuncAtt(fitCorr, colorTot, styleTot, lineSize);
    SetFuncAtt(fitCorrSig, colorSig, styleSig, lineSize);
    SetFuncAtt(fitCorrBg, colorBg, styleBg, lineSize);

    gSystem->mkdir(Form("%s/fits/corr/single/",sPath.Data()),1);

    gCanSingle->cd();
    histCorr->DrawCopy();
    fitCorr->Draw("same");
    fitCorrSig->Draw("same");
    fitCorrBg->Draw("same");
    gCanSingle->SaveAs(Form("%s/fits/corr/single/histCorr_%s.%s",sPath.Data(),sListName.Data(),gOutFormat.Data()),gOutFormat.Data());
    gCanSingle->Clear();

    if(padIndex > 0) {
        gCanOverCorr->cd(padIndex);
        SetPad();
        SetHistOff(histCorr);
        histCorr->DrawCopy();
        fitCorr->Draw("same");
        fitCorrSig->Draw("same");
        fitCorrBg->Draw("same");
    }

    // ========== Processing histFrac =======================================

    TH1D* histMass_fracSig = (TH1D*) list->FindObject("histMass_fracSig");
    if(!histMass_fracSig) { printf("ERROR: Object '%s' not found!\n", "histMass_fracSig"); list->ls(); return kFALSE; }

    TH1D* histMass_fracBg = (TH1D*) list->FindObject("histMass_fracBg");
    if(!histMass_fracBg) { printf("ERROR: Object '%s' not found!\n", "histMass_fracBg"); list->ls(); return kFALSE; }

    TF1* fitMass_fracSig = (TF1*) list->FindObject("fitMass_fracSig");
    if(!fitMass_fracSig) { printf("ERROR: Object '%s' not found!\n", "fitMass_fracSig"); list->ls(); return kFALSE; }

    TF1* fitMass_fracBg = (TF1*) list->FindObject("fitMass_fracBg");
    if(!fitMass_fracBg) { printf("ERROR: Object '%s' not found!\n", "fitMass_fracBg"); list->ls(); return kFALSE; }

    // Double_t dXlow = fitCorr->GetXmin();
    // Double_t dXhigh = fitCorr->GetXmax();

    histMass_fracSig->SetTitle(sTitle.Data());
    histMass_fracSig->GetXaxis()->SetRangeUser(dXlow,dXhigh);
    histMass_fracSig->SetStats(0);
    histMass_fracSig->SetMinimum(-0.05);
    histMass_fracSig->SetMaximum(1.05);

    SetHistAtt(histMass_fracSig, Color_t(colorSig+1), kFullCircle, markSize);
    SetHistAtt(histMass_fracBg, Color_t(colorBg+2), kFullCircle, markSize);

    SetFuncAtt(fitMass_fracSig, colorSig, kSolid, lineSize);
    SetFuncAtt(fitMass_fracBg, colorBg, kSolid, lineSize);

    gSystem->mkdir(Form("%s/fits/mass/single/",sPath.Data()),1);

    gCanSingle->cd();
    histMass_fracSig->DrawCopy();
    histMass_fracBg->DrawCopy("same");
    fitMass_fracSig->Draw("same");
    fitMass_fracBg->Draw("same");
    gCanSingle->SaveAs(Form("%s/fits/mass/single/mass_%s.%s",sPath.Data(),sListName.Data(),gOutFormat.Data()),gOutFormat.Data());
    gCanSingle->Clear();

    if(padIndex > 0) {
        gCanOverMass->cd(padIndex);
        SetPad();
        SetHistOff(histMass_fracSig);
        histMass_fracSig->DrawCopy();
        histMass_fracBg->DrawCopy("same");
        fitMass_fracSig->Draw("same");
        fitMass_fracBg->Draw("same");
    }


    // ========== Processing histFrac =======================================

    return kTRUE;
}

void PrepareCanvas(Int_t iNumPt = 8)
{
    gCanSingle = new TCanvas("gCanSingle","gCanSingle", 400,400);
    // gCanSingle->cd(1);

    Int_t iCols = 3;
    Int_t iRow = TMath::Ceil(Double_t(iNumPt)/iCols);

    gCanOverMass = new TCanvas("gCanOverMass","gCanOverMass", 700,700);
    gCanOverMass->Divide(iCols,iRow);

    gCanOverCorr = new TCanvas("gCanOverCorr","gCanOverCorr", 700,700);
    gCanOverCorr->Divide(iCols,iRow);
}

void SetPad()
{
    gPad->SetBottomMargin(0.2);
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.13);
}

void SetHistOff(TH1* hist)
{
    hist->GetXaxis()->SetNdivisions(505,1);

    hist->SetLabelFont(43,"XY");
    hist->SetLabelSize(12,"XY");

    hist->SetTitleFont(43,"XY");
    hist->SetTitleSize(12,"XY");

    hist->SetTitleOffset(4.3,"X");
    // hist->SetTitleOffset(4.3,"t");
    // hist->SetTitleFont(43,"t");
    // hist->SetTitleSize(15,"t");

    // hist->SetTitleOffset(2.2,"XY");
}

void SetFuncAtt(TF1* func, Color_t color, Style_t style, Int_t width)
{
    func->SetLineColor(color);
    func->SetLineStyle(style);
    func->SetLineWidth(width);
}

void SetHistAtt(TH1* hist, Color_t color, Style_t markStyle, Double_t markSize)
{
    hist->SetMarkerColor(color);
    hist->SetMarkerStyle(markStyle);
    hist->SetMarkerSize(markSize);

    hist->SetLineColor(color);
    // hist->SetLineStyle(lineStyle);
    // hist->SetLineWidth(lineWidth);
}
