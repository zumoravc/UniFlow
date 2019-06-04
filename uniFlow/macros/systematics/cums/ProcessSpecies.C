Int_t iNumCent = 6;
std::vector<TString> vecSyst = {};

// Int_t iNumCent = 6;
// TString path = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/PbPb/syst_6815/" + sSpecies + "/";
// TString baseline = "gap00";
// std::vector<TString> vecSyst = {"CL1"};
//

void ProcessSpecies(
    std::vector<TString> vecSpecies,
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results//PbPb/cums/syst_6815/",
    TString sPathOut = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results//PbPb/cums/syst/",
    TString sBaseline = "gap00"
);

void ProcessSpeciesAll(
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_6815/",
    TString sPathOut = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_fluct/",
    TString sBaseline = "gap00")
{
    ProcessSpecies({"Charged","Pion","Kaon","Proton","K0s","Lambda"},sPath,sPathOut,sBaseline);
}


void ProcessSpecies(
    TString sSpecies = "Charged",
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_6815/",
    TString sPathOut = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst/",
    TString sBaseline = "gap00"
)
{
    TString path = sPath + sSpecies + "/";
    TString pathOut = sPathOut + sSpecies + "/";

    for(Int_t iC(0); iC < iNumCent; ++iC) {

        TString hist = Form("%s_hFlow4_harm2_gap-10_cent%d",sSpecies.Data(),iC);

        TMacro single = TMacro(gSystem->ExpandPathName("~/Codes/ALICE/Flow/uniFlow/macros/systematics/cums/ProcessSingle.C"));

        Double_t dFitXmin = 3.0;
        Double_t dFitXminPID3sigma = 1.0;

        if(sSpecies.EqualTo("Charged")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"FB768",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),dFitXmin));
        } else if(sSpecies.EqualTo("Pion") || sSpecies.EqualTo("Kaon") || sSpecies.EqualTo("Proton")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"FB768",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma",sBaseline.Data(),dFitXminPID3sigma));
        } else if(sSpecies.EqualTo("K0s") || sSpecies.EqualTo("Lambda")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"V0sCPA099",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"V0sCrossFind1",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"V0sDaugDCA3",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"V0sDaugPt02",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"V0sDecRad1",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"V0sDecRad10",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"V0sFinderOn",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"V0sPVDCA3",sBaseline.Data(),-1.0));
        } else if(sSpecies.EqualTo("Phi")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"FB768",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma",sBaseline.Data(),-1.0));
        } else {
            return kFALSE;
        }
    }
}

void ProcessSpecies(
    std::vector<TString> vecSpecies,
    TString sPath,
    TString sPathOut,
    TString sBaseline
)
{
    Int_t iNumSp = vecSpecies.size();
    for(Int_t iS(0); iS < iNumSp; ++iS) {
        ProcessSpecies(vecSpecies.at(iS),sPath,sPathOut,sBaseline);
    }
}
