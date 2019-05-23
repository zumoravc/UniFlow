// Macro for transfering TH1* to TGraphErrors object
#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/utils/Utils.cxx>
#include "TH1.h"
#include "TGraphErrors.h"

void SetXerrs(TGraphErrors* graph, Double_t dErrX); // Set graphs X-axis errors to Ex value

TGraphErrors* HistoToGraph(TH1* histo, Double_t dErrX = -1)
{
    if(!histo) {
        Utils::Error("Input histo not found!","HistoToGraph");
        return nullptr;
    }

    TGraphErrors* err = new TGraphErrors(histo);
    SetXerrs(err,dErrX);
    return err;
}

void SetXerrs(TGraphErrors* graph, Double_t dErrX)
{
    if(!graph) {
        Utils::Error("Input graph not found!","SetXerrs");
        return;
    }

    if(dErrX < 0.0) {
        Utils::Debug("Value of ErrX < 0! Skipping!",1,"SetXerrs");
        return;
    }

    Int_t np = graph->GetN();
    for (Int_t i(0); i < np; ++i) {
        Double_t dErrY = graph->GetErrorY(i);
        graph->SetPointError(i, dErrX, dErrY);
    }
}
