#ifndef UTILS_CXX
#define UTILS_CXX

#include <TString.h>
#include "Utils.h"

Int_t Utils::fgDebugLvl = 0;

void Utils::Message(EMsg type, TString sMsg, TString sMethod)
{
    const TString Msg_Header[] = { "Debug", "Info", "Warning", "Error", "Fatal"};
    const TString Msg_Color[] = {"95m", "96m", "93m", "91m", "91m"};

    TString sHeader = Msg_Header[type];
    if(sMethod.Length()) { sHeader += "-"; sHeader += sMethod; }
    printf("\033[%s%s : %s\033[0m\n", Msg_Color[type].Data(), sHeader.Data(), sMsg.Data());
};

#endif
