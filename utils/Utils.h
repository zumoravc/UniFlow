#ifndef UTILS_H
#define UTILS_H

#include "TString.h"
class TH1D;
class TLegend;

class Utils
{
    public:
        enum PosLegend {kLegTopLeft = 1, kLegTopRight, kLegBotLeft, kLegBotRight};

        static TH1D* DivideHistos(TH1D* nom, TH1D* denom, Bool_t bCor = kFALSE);

        static void SetDebugLevel(Int_t iLvl = 1) { fgDebugLvl = iLvl; }

        static void Debug(TString sMsg, Int_t iLvl = 0, TString sMethod = "") { if(fgDebugLvl >= iLvl) Message(EMsg(kDebug), Form("%s (lvl %d)", sMsg.Data(), iLvl), sMethod); }
        static void Info(TString sMsg, TString sMethod = "") { Message(EMsg(kInfo), sMsg, sMethod); }
        static void Warning(TString sMsg, TString sMethod = "") { Message(EMsg(kWarning), sMsg, sMethod); }
        static void Error(TString sMsg, TString sMethod = "") { Message(EMsg(kError), sMsg, sMethod); }
        static void Fatal(TString sMsg, TString sMethod = "") { Message(EMsg(kFatal), sMsg, sMethod); }

        static TLegend* MakeLegend(PosLegend pos);

    private:
        enum EMsg { kDebug = 0, kInfo, kWarning, kError, kFatal, kNum};
        static Int_t fgDebugLvl;
        static void Message(EMsg type, TString sMsg, TString sMethod = "");

};
#endif
