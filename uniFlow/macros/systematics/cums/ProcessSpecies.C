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
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results//PbPb/cums/syst_6815/",
    TString sBaseline = "gap00"
);

void ProcessSpeciesAll(
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_6815/",
    TString sBaseline = "gap00")
{
    ProcessSpecies({"Charged","Pion","Kaon","Proton","K0s","Lambda","Phi"},sPath,sBaseline);
}


void ProcessSpecies(
    TString sSpecies = "Charged",
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_6815/",
    TString sBaseline = "gap00"
)
{
    TString path = sPath + sSpecies + "/";

    for(Int_t iC(0); iC < iNumCent; ++iC) {

        TString hist = Form("%s_hFlow4_harm2_gap-10_cent%d",sSpecies.Data(),iC);

        TMacro single = TMacro(gSystem->ExpandPathName("~/Codes/ALICE/Flow/uniFlow/macros/systematics/cums/ProcessSingle.C"));

        Double_t dFitXmin = 3.0;
        Double_t dFitXminPID3sigma = 1.0;

        if(sSpecies.EqualTo("Charged")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"CL1",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"FB768",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"PVz8",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"TPCcls90",sBaseline.Data(),dFitXmin));
        } else if(sSpecies.EqualTo("Pion") || sSpecies.EqualTo("Kaon") || sSpecies.EqualTo("Proton")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"CL1",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"FB768",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"PVz8",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"TPCcls90",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"PID3sigma",sBaseline.Data(),dFitXminPID3sigma));
        } else if(sSpecies.EqualTo("K0s") || sSpecies.EqualTo("Lambda")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"CL1",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"PVz8",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"TPCcls90",sBaseline.Data(),dFitXmin));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"V0sCPA099",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"V0sCrossFind1",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"V0sDaugDCA3",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"V0sDaugPt02",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"V0sDecRad1",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"V0sDecRad10",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"V0sFinderOn",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"V0sPVDCA3",sBaseline.Data(),-1.0));
        } else if(sSpecies.EqualTo("Phi")) {
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"CL1",sBaseline.Data(),-1.0));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"FB768",sBaseline.Data(),-1.0.q));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"PVz8",sBaseline.Data(),-1.0.q));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"TPCcls90",sBaseline.Data(),-1.0.q));
            single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),"PID3sigma",sBaseline.Data(),-1.0.q));
        } else {
            return kFALSE;
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

// Bool_t SetArgs(TString sp)
// {
//     if(sSpecies.EqualTo("Charged")) {
//         vecSyst = {"CL1","FB768","PVz8","TPCcls90"};
//     } else if(sSpecies.EqualTo("Pion")) {
//         vecSyst = {"CL1","FB768","PID3sigma","PVz8","TPCcls90"};
//     } else if(sSpecies.EqualTo("Kaon")) {
//         vecSyst = {"CL1","FB768","PID3sigma","PVz8","TPCcls90"};
//     } else if(sSpecies.EqualTo("Proton")) {
//         vecSyst = {"CL1","FB768","PID3sigma","PVz8","TPCcls90"};
//     } else if(sSpecies.EqualTo("K0s")) {
//         vecSyst = {
//             "CL1","PVz8","TPCcls90","V0sCPA099","V0sCrossFind1","V0sDaugDCA3","V0sDaugPt02",
//             "V0sDecRad1","V0sDecRad10","V0sFinderOn","V0sPVDCA3"
//         };
//     } else if(sSpecies.EqualTo("Lambda")) {
//         vecSyst = {
//             "CL1","PVz8","TPCcls90","V0sCPA099","V0sCrossFind1","V0sDaugDCA3","V0sDaugPt02",
//             "V0sDecRad1","V0sDecRad10","V0sFinderOn","V0sPVDCA3"
//         };
//     } else if(sSpecies.EqualTo("Phi")) {
//         // iNumCent = 5;
//         vecSyst = {"CL1","FB768","PID3sigma","PVz8","TPCcls90"};
//     } else {
//         return kFALSE;
//     }
//
//     return kTRUE;
// }
