#ifndef UTILS_CXX
#define UTILS_CXX

#include "Utils.h"
#include <TString.h>
#include <TLegend.h>

Int_t Utils::fgDebugLvl = 0;

void Utils::Message(EMsg type, TString sMsg, TString sMethod)
{
    const TString Msg_Header[] = { "Debug", "Info", "Warning", "Error", "Fatal"};
    const TString Msg_Color[] = {"95m", "96m", "93m", "91m", "91m"};

    TString sHeader = Msg_Header[type];
    if(sMethod.Length()) { sHeader += "-"; sHeader += sMethod; }
    printf("\033[%s%s : %s\033[0m\n", Msg_Color[type].Data(), sHeader.Data(), sMsg.Data());
};

// ===========================================================================
TLegend* Utils::MakeLegend(PosLegend pos)
{
    TLegend* leg = nullptr;

    switch(pos) {
        case kLegTopLeft: leg = new TLegend(0.14,0.65,0.5,0.88); break;
        case kLegTopRight: leg = new TLegend(0.58,0.65,0.88,0.88); break;
        case kLegBotLeft: leg = new TLegend(0.14,0.12,0.5,0.4); break;
        case kLegBotRight: leg = new TLegend(0.58,0.12,0.88,0.4); break;
        default: return nullptr;
    }

    leg->SetBorderSize(0);
    leg->SetFillColorAlpha(0,0.0);

    return leg;
}
// ===========================================================================
TH1D* Utils::DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor)
{
    if(!nom || !denom) { printf("ERR: either of the histos does not exists\n"); return 0x0; }

    Int_t binsNom = nom->GetNbinsX();
    Int_t binsDenom = denom->GetNbinsX();

    if(binsNom != binsDenom) { Utils::Warning("Different # of bins\n"); }
    Int_t binsMin = (binsNom < binsDenom ? binsNom : binsDenom);

    TH1D* ratio = (TH1D*) nom->Clone(Form("Ratio_%s_%s",nom->GetName(),denom->GetName()));
    ratio->Reset();

    for(Int_t iBin(1); iBin < binsDenom+1; ++iBin) {
        if(iBin > binsMin) break;

        Double_t dContNom = nom->GetBinContent(iBin);
        Double_t dErrNom = nom->GetBinError(iBin);
        Double_t dContDenom = denom->GetBinContent(iBin);
        Double_t dErrDenom = denom->GetBinError(iBin);

        if(dContDenom == 0.0) { continue; }

        Double_t dContRatio = dContNom / dContDenom;

        Double_t dContribNom = dErrNom / dContDenom;
        Double_t dContribDenom = -1.0 * dContNom * dErrDenom * TMath::Power(dContDenom,-2.0);

        Double_t dErrRatioSq = 0.0;
        if(bCor) {
            dErrRatioSq = TMath::Power(dContribNom + dContribDenom,2.0);
        } else {
            dErrRatioSq = TMath::Power(dContribNom,2.0) + TMath::Power(dContribDenom,2.0);
        }

        ratio->SetBinContent(iBin,dContRatio);
        ratio->SetBinError(iBin,TMath::Sqrt(dErrRatioSq));
    }

    return ratio;
}

#endif
