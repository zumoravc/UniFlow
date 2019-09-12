// Running macro for processing systematics
// Setup for additional PID syst. checks

Int_t iNumCent = 6;
std::vector<TString> vecSyst = {};

void ProcessSpecies(
    std::vector<TString> vecSpecies,
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results//PbPb/cums/syst_6883/",
    TString sPathOut = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results//PbPb/cums/syst/",
    TString sBaseline = "gap00"
);

void ProcessSpeciesAll(
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_6883/",
    TString sPathOut = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results//PbPb/cums/syst/",
    TString sBaseline = "gap00")
{
    ProcessSpecies({"Pion","Kaon","Proton","Phi"},sPath,sPathOut,sBaseline);
}


void ProcessSpecies(
    TString sSpecies = "Pion",
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_6883/",
    TString sPathOut= "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst/",
    TString sBaseline = "gap00"
)
{
    TString path = sPath + sSpecies + "/";
    TString pathOut = sPathOut + sSpecies + "/";

    for(Int_t iC(0); iC < iNumCent; ++iC) {

        TString hist = Form("%s_hFlow4_harm2_gap-10_cent%d",sSpecies.Data(),iC);

        TMacro single = TMacro(gSystem->ExpandPathName(Form("%s/macros/systematics/cums/ProcessSingle.C",gSystem->Getenv("UNIFLOW"))));

        Double_t dFitXminPID3sigma = 1.0;

        if(sSpecies.EqualTo("Pion")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"Bayes90",sBaseline.Data(),1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID2sigma",sBaseline.Data(),1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma",sBaseline.Data(),1.0));
        }

        if(sSpecies.EqualTo("Kaon")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"Bayes90",sBaseline.Data(),1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID2sigma",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma",sBaseline.Data(),-1.0 ));
        }

        if(sSpecies.EqualTo("Proton")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"Bayes90",sBaseline.Data(),2.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID2sigma",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma",sBaseline.Data(),-1.0));

            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID2sigma_anti",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma_anti",sBaseline.Data(),-1.0));
        }

        if(sSpecies.EqualTo("Phi")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"Bayes90",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID2sigma",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma",sBaseline.Data(),-1.0));
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
