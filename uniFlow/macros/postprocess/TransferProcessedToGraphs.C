// Macro for transfering TH1* to TGraphErrors object

#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/macros/postprocess/HistoToGraph.C>
#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/macros/utils/Utils.cxx>

#include <vector>
#include "TFile.h"
#include <TString.h>

void TransferProcessedToGraphs(
    const std::vector<TString> vecHistos = {"hFlow2_harm2_gap-10","hFlow2_harm2_gap08", "hFlow4_harm2_gap-10"},
    const std::vector<TString> vecSpecies = {"Charged","Pion","Kaon","K0s","Proton","Lambda", "Phi"},
    const Int_t iNumCent = 1,
    const TString sInputFile = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/pPb/2464/V0A_020/gap_all/Processed.root",
    const TString sOutputFile = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/pPb/2464/V0A_020/gap_all/Graphs.root",
    const TString sFileOutMode = "RECREATE",
    const Bool_t bIsSyst = kFALSE
)
{
    TFile* fileIn = TFile::Open(sInputFile.Data(),"READ");
    if(!fileIn) {
        Utils::Error("file input not open!");
        return;
    }

    TFile* fileOut = TFile::Open(sOutputFile.Data(),sFileOutMode.Data());
    if(!fileOut) {
        Utils::Error("file output not open!");
        fileIn->Close();
        return;
    }

    const Int_t iNumHistos = vecHistos.size();
    const Int_t iNumSpecies = vecSpecies.size();

    for(Int_t iH(0); iH < iNumHistos; ++iH) {

        for(Int_t iS(0); iS < iNumSpecies; ++iS) {

            for(Int_t iC(0); iC < iNumCent; ++iC) {

                TString sHistoName = Form("%s_%s_cent%d",vecSpecies.at(iS).Data(), vecHistos.at(iH).Data(), iC);

                TH1* hist = (TH1*) fileIn->Get(sHistoName.Data());
                if(!hist) {
                    Utils::Error(Form("histo '%s' not found!", sHistoName.Data()));
                    fileIn->ls();
                    fileIn->Close();
                    fileOut->Close();
                    return;
                }

                TGraphErrors* graph = HistoToGraph(hist);
                if(!graph) {
                    Utils::Error("Graph transfer failed!");
                } else {
                    TString sHistoOutName = sHistoName;
                    if(bIsSyst) { sHistoOutName += "_syst"; }

                    graph->SetName(sHistoOutName.Data());
                    fileOut->cd();
                    graph->Write();
                }
            }
        }
    }

    fileIn->Close();
    fileOut->Close();

    return;
}

void TransferProcessedToGraphsWithSyst(
    const std::vector<TString> vecHistos = {"hFlow2_harm2_gap-10", "hFlow2_harm2_gap08", "hFlow4_harm2_gap-10"},
    const std::vector<TString> vecSpecies = {"Charged","Pion","Kaon","K0s","Proton","Lambda", "Phi"},
    const Int_t iNumCent = 6,
    const TString sInputFile = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/pPb/2464/V0A_020/gap_all/Processed.root",
    const TString sInputFileSyst = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/pPb/syst_2517/"
    const TString sOutputFile = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/pPb/2464/V0A_020/gap_all/Cums_pPb.root",
)
{
    TransferProcessedToGraphs(vecHistos, vecSpecies, iNumCent, sInputFile, sOutputFile, "RECREATE",0);
    TransferProcessedToGraphs(vecHistos, vecSpecies, iNumCent, sInputFileSyst, sOutputFile, "UPDATE",1);
};
