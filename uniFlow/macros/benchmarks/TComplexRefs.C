#include <iostream>
#include "TComplex.h"
#include "TStopwatch.h"
#include "TH1.h"

TComplex fComp = TComplex(1.0,-2.0,0);

TComplex CompCopy() { return fComp; }
TComplex& CompRef() { return fComp; }

void TComplexRefs(ULong64_t iNumTot = 1e9, Int_t iNumRuns = 1)
{
    // Int_t iNumRuns = 1;
    // ULong64_t iNumTot = 1e9;

    TStopwatch sw = TStopwatch();

    TH1D* histCopy = new TH1D("histCopy","histCopy",50,0,5);
    TH1D* histRef = new TH1D("histRef","histRef",50,0,5);

    for(Int_t iRun(0); iRun < iNumRuns; ++iRun) {

        // returning copy
        sw.Start(kTRUE);
        for(ULong64_t iNum(0); iNum < iNumTot; ++iNum) {
            TComplex comp = CompCopy();
            // std::cout << comp << std::endl;

        }
        sw.Stop();

        histCopy->Fill(sw.CpuTime());
        // histCopy->Fill(sw.RealTime());

        std::cout << "### Copy ### Run " << iRun << " | CPU : "<< sw.CpuTime() << " | Real : " << sw.RealTime() << " (after " << iNumTot << " calls) " << std::endl;

        // returning reference
        sw.Start(kTRUE);
        for(ULong64_t iNum(0); iNum < iNumTot; ++iNum) {
            TComplex comp = CompRef();
            // std::cout << comp << std::endl;
        }
        sw.Stop();

        histRef->Fill(sw.CpuTime());
        // histRef->Fill(sw.RealTime());

        std::cout << "### Ref ### Run " << iRun << " | CPU : "<< sw.CpuTime() << " | Real : " << sw.RealTime() << " (after " << iNumTot << " calls) " << std::endl;
    }

    std::cout << "Average (CPU time) : Copy : " << histCopy->GetMean() << " ± " << histCopy->GetRMS() << " | Refs : " << histRef->GetMean() << " ± "<< histRef->GetRMS() << " (ratio " << histCopy->GetMean() / histRef->GetMean() << ")"<< std::endl;

    histCopy->SetLineColor(kRed);
    histCopy->SetMarkerColor(kRed);

    histRef->SetLineColor(kBlack);
    histRef->SetMarkerColor(kBlack);

    histCopy->Draw("same");
    histRef->Draw("same");
}
