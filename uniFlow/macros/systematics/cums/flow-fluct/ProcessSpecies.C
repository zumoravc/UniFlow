Int_t iNumCent = 4;
std::vector<TString> vecSyst = {};

// Int_t iNumCent = 6;
// TString path = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/PbPb/syst_6815/" + sSpecies + "/";
// TString baseline = "gap00";
// std::vector<TString> vecSyst = {"CL1"};
//

void ProcessSpecies(
    std::vector<TString> vecSpecies,
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_fluct/",
    TString sPathOut = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_fluct/",
    TString sBaseline = "gap00"
);

void ProcessSpeciesAll(
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_fluct/",
    TString sPathOut = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_fluct/",
    TString sBaseline = "gap00")
{
    ProcessSpecies({"Charged","Pion","Kaon","Proton","K0s","Lambda"},sPath,sPathOut,sBaseline);
}


void ProcessSpecies(
    TString sSpecies = "Charged",
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_fluct/",
    TString sPathOut = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_fluct/",
    TString sBaseline = "gap00"
)
{
    TString path = sPath + sSpecies + "/";
    TString pathOut = sPathOut + sSpecies + "/";

    for(Int_t iC(0); iC < iNumCent; ++iC) {

        TString hist = Form("%s_rel_cent%d",sSpecies.Data(),iC);

        TMacro single = TMacro(gSystem->ExpandPathName("~/Codes/ALICE/Flow/uniFlow/macros/systematics/cums/flow-fluct/ProcessSingle.C"));

        Double_t dFitXmin = 3.0;
        Double_t dFitXminPID3sigma = 1.0;

        Bool_t bCor = 1;

        if(sSpecies.EqualTo("Charged")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),2.5, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"FB768",sBaseline.Data(),-1.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),-1.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),-1.0, bCor));
        } else if(sSpecies.EqualTo("Pion")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"FB768",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PID2sigma",sBaseline.Data(),0.0, bCor));
            // single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"Bayes90",sBaseline.Data(),dFitXminPID3sigma, bCor));
        } else if(sSpecies.EqualTo("Kaon")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"FB768",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PID2sigma",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"Bayes90",sBaseline.Data(),0.0, bCor));
        } else if(sSpecies.EqualTo("Proton")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"FB768",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"Bayes90",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PID2sigma_anti",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma_anti",sBaseline.Data(),0.0, bCor));
        } else if(sSpecies.EqualTo("K0s") || sSpecies.EqualTo("Lambda")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"V0sCPA099",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"V0sCrossFind1",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"V0sDaugDCA3",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"V0sDaugPt02",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"V0sDecRad1",sBaseline.Data(),0.0, bCor));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"V0sDecRad10",sBaseline.Data(),0.0, bCor));
            // single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"V0sFinderOn",sBaseline.Data(),-1.0, bCor));
            // single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"V0sPVDCA3",sBaseline.Data(),-1.0, bCor));
        // } else if(sSpecies.EqualTo("Phi")) {
            // single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),-1.0, bCor));
            // single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"FB768",sBaseline.Data(),-1.0, bCor));
            // single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),-1.0, bCor));
            // single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),-1.0, bCor));
            // single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f,%d",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma",sBaseline.Data(),-1.0, bCor));
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
