// Macro for comparing p-Pb results of v2{4} with v2{2}^{sub} (Preliminaries)

#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/macros/utils/Utils.cxx>
#include <TString.h>
#include <TFile.h>

void PlotPPb(
    TString sSpecies = "Pion",
    Int_t iCent = 0
)
{
    TString sPathOut = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/pPb/2464/V0A_020/gap_all/plots_compCumsSub/";

    TString sFileSub = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/v22_subt_pPb_QM/v2-syst-final.root";
    TString sFileCums = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/pPb/2464/V0A_020/gap_all/Graphs.root";

    TString sHistoSub = Form("graphPoints_%s_cent%d", sSpecies.Data(), iCent);
    TString sHistoCumTwo = Form("%s_hFlow2_harm2_gap-10_cent%d", sSpecies.Data(), iCent);
    TString sHistoCumFour = Form("%s_hFlow4_harm2_gap-10_cent%d", sSpecies.Data(), iCent);

    Color_t col[] = { kBlack, kRed, kBlue};
    Int_t markers[] = { kOpenCircle, kOpenCircle, kOpenCircle};

    /////////

    TFile* fileSub = TFile::Open(sFileSub.Data(), "READ");
    if(!fileSub) {  Utils::Error("fileSub not found!"); return; }
    TGraphErrors* grSub = (TGraphErrors*) fileSub->Get(sHistoSub.Data());
    if(!grSub) { Utils::Error(Form("grSub '%s' not found!",sFileSub.Data())); fileSub->ls(); return; }
    grSub->SetMarkerStyle(markers[0]);
    grSub->SetMarkerColor(col[0]);
    grSub->SetLineColor(col[0]);

    TFile* fileCums = TFile::Open(sFileCums.Data(), "READ");
    if(!fileCums) {  Utils::Error("fileCums not found!"); return; }

    TGraphErrors* grCumTwo = (TGraphErrors*) fileCums->Get(sHistoCumTwo);
    if(!grCumTwo) { Utils::Error(Form("grCumTwo '%s' not found!",sHistoCumTwo.Data())); fileCums->ls(); return; }
    grCumTwo->SetMarkerStyle(markers[1]);
    grCumTwo->SetMarkerColor(col[1]);
    grCumTwo->SetLineColor(col[1]);

    TGraphErrors* grCumFour = (TGraphErrors*) fileCums->Get(sHistoCumFour);
    if(!grCumFour) { Utils::Error(Form("grCumFour '%s' not found!",sHistoCumFour.Data())); fileCums->ls(); return; }
    grCumFour->SetMarkerStyle(markers[2]);
    grCumFour->SetMarkerColor(col[2]);
    grCumFour->SetLineColor(col[2]);

    TCanvas* can = new TCanvas("can","can", 600,600);
    can->cd();

    TLegend* leg = Utils::MakeLegend(Utils::kLegTopLeft);
    leg->SetHeader(Form("%s (stat. only)",sSpecies.Data()));
    leg->AddEntry(grCumTwo, "v_{2}{2}", "pl");
    leg->AddEntry(grSub, "v_{2}^{sub}{2,|#Delta#eta|>0.4}", "pl");
    leg->AddEntry(grCumFour, "v_{2}{4}", "pl");

    TH1* frame = (TH1*) gPad->DrawFrame(0.0,0.0,6.0,0.5);

    grSub->Draw("same pe");
    grCumTwo->Draw("same pe");
    grCumFour->Draw("same pe");
    leg->Draw();

    gSystem->mkdir(sPathOut.Data(),1);
    can->SaveAs(Form("%s/%s_cent%d.pdf",sPathOut.Data(), sSpecies.Data(), iCent), "pdf");

    fileSub->Close();
    fileCums->Close();
}

void PlotAll(Int_t iCent = 0)
{
    PlotPPb("Charged",iCent);
    PlotPPb("Pion",iCent);
    PlotPPb("Kaon",iCent);
    PlotPPb("Proton",iCent);
    PlotPPb("K0s",iCent);
    PlotPPb("Lambda",iCent);
    PlotPPb("Phi",iCent);
}
