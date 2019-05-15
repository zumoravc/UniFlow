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


#endif
