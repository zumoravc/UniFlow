Int_t iNumCent = 6;
std::vector<TString> vecSyst = {};

// Int_t iNumCent = 6;
// TString path = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/PbPb/syst_6815/" + sSpecies + "/";
// TString baseline = "gap00";
// std::vector<TString> vecSyst = {"CL1"};
//

Bool_t SetArgs(TString sp);

void ProcessSpecies(
    std::vector<TString> vecSpecies,
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/PbPb/syst_6815/",
    TString sBaseline = "gap00"
);

void ProcessSpeciesAll(
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/PbPb/syst_6815/",
    TString sBaseline = "gap00")
{
    ProcessSpecies({"Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"},sPath,sBaseline);
}


void ProcessSpecies(
    TString sSpecies = "Charged",
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/cums/PbPb/syst_6815/",
    TString sBaseline = "gap00"
)
{
    if(!SetArgs(sSpecies)) {
        printf("ERROR: invalid species!\n");
        return;
    }

    TString path = sPath + sSpecies + "/";

    Int_t iNumSyst = vecSyst.size();

    for(Int_t iS(0); iS < iNumSyst; ++iS) {

        TString syst = vecSyst.at(iS);

        for(Int_t iC(0); iC < iNumCent; ++iC) {

            TString hist = Form("%s_hFlow4_harm2_gap-10_cent%d",sSpecies.Data(),iC);

            TMacro single = TMacro(gSystem->ExpandPathName("~/Codes/ALICE/Flow/uniFlow/macros/syst_cums/ProcessSingle.C"));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\"",hist.Data(),path.Data(),syst.Data(),sBaseline.Data()));
        }
    }
}

void ProcessSpecies(
    std::vector<TString> vecSpecies,
    TString sPath,
    TString sBaseline
)
{
    Int_t iNumSp = vecSpecies.size();
    for(Int_t iS(0); iS < iNumSp; ++iS) {
        ProcessSpecies(vecSpecies.at(iS),sPath,sBaseline);
    }
}

Bool_t SetArgs(TString sp)
{
    if(sp.EqualTo("Charged")) {
        vecSyst = {"CL1","FB768","PVz8","TPCcls90"};
    } else if(sp.EqualTo("Pion")) {
        vecSyst = {"CL1","FB768","PID3sigma","PVz8","TPCcls90"};
    } else if(sp.EqualTo("Kaon")) {
        vecSyst = {"CL1","FB768","PID3sigma","PVz8","TPCcls90"};
    } else if(sp.EqualTo("Proton")) {
        vecSyst = {"CL1","FB768","PID3sigma","PVz8","TPCcls90"};
    } else if(sp.EqualTo("K0s")) {
        vecSyst = {
            "CL1","PVz8","TPCcls90","V0sCPA099","V0sCrossFind1","V0sDaugDCA3","V0sDaugPt02",
            "V0sDecRad1","V0sDecRad10","V0sFinderOn","V0sPVDCA3"
        };
    } else if(sp.EqualTo("Lambda")) {
        vecSyst = {
            "CL1","PVz8","TPCcls90","V0sCPA099","V0sCrossFind1","V0sDaugDCA3","V0sDaugPt02",
            "V0sDecRad1","V0sDecRad10","V0sFinderOn","V0sPVDCA3"
        };
    } else if(sp.EqualTo("Phi")) {
        // iNumCent = 5;
        vecSyst = {"CL1","FB768","PID3sigma","PVz8","TPCcls90"};
    } else {
        return kFALSE;
    }

    return kTRUE;
}
