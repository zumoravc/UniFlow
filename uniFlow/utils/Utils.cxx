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

  // if(binsNom != binsDenom) { printf("ERR: Different # of bins\n"); return 0x0; }

  TH1D* ratio = (TH1D*) nom->Clone(Form("Ratio_%s_%s",nom->GetName(),denom->GetName()));
  ratio->Reset();

  Double_t dContNom = 0, dErrNom = 0;
  Double_t dContDenom = 0, dErrDenom = 0;
  Double_t dContRatio = 0, dErrRatio = 0;
  for(Short_t iBin(1); iBin < binsDenom+1; iBin++)
  {
    if(iBin > binsNom) break;

    dContNom = nom->GetBinContent(iBin);
    dErrNom = nom->GetBinError(iBin);
    dContDenom = denom->GetBinContent(iBin);
    dErrDenom = denom->GetBinError(iBin);

    dContRatio =  dContNom / dContDenom;
    dErrRatio = TMath::Power(dErrNom/dContDenom, 2) + TMath::Power( dErrDenom*dContNom/(dContDenom*dContDenom), 2);
    // printf("Err (before) : %g | ", TMath::Sqrt(dErrRatio));

    if(bCor) dErrRatio -= (2*dContNom*dErrDenom*dErrNom/TMath::Power(dContDenom,3));
    // printf("(after) : %g\n", TMath::Sqrt(dErrRatio));

    ratio->SetBinContent(iBin,dContRatio);
    ratio->SetBinError(iBin,TMath::Sqrt(dErrRatio));
  }

  return ratio;
}

#endif
