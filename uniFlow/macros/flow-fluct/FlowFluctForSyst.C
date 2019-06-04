#include </Users/vpacik/Codes/ALICE/Flow/uniFlow/utils/Utils.cxx>

void FlowFluctForSyst(
    TString sSpecies = "Charged",
    TString sPath = "/Users/vpacik/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/syst_fluct",
    Int_t iNumCent = 4,
    Bool_t bCor = 1,
    TString sHistoTwo = "hFlow2_harm2_gap-10",
    TString sHistoFour = "hFlow4_harm2_gap-10",
    TString sOutFileName = "FlowFluct"
)
{
    // FlowFluct(TString sSpecies = "Charged", Int_t iCent = 0, TString sPath = "~/Codes/ALICE/Flow/uniFlow/results/PbPb/cums/MC_2210/species_direct/", Bool_t bCor = 0);
    TMacro FlowFluct = TMacro(gSystem->ExpandPathName("~/Codes/ALICE/Flow/uniFlow/macros/flow-fluct/FlowFluct.C"));
    // FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, sPath, bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));

    for(Int_t iCent = 0; iCent < iNumCent; ++iCent) {

        // Int_t iCent = 0;

        if(sSpecies.EqualTo("Charged")) {
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"gap00"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"CL1"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"FB768"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"PVz8"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"TPCcls90"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
        } else if(sSpecies.EqualTo("Pion") || sSpecies.EqualTo("Kaon") || sSpecies.EqualTo("Proton")) {
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"gap00"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"CL1"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"FB768"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"PVz8"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"TPCcls90"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"Bayes90"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"PID2sigma"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"PID3sigma"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));

            if(sSpecies.EqualTo("Proton")) {
                FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"PID3sigma_anti"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
                FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"PID2sigma_anti"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            }

        } else if(sSpecies.EqualTo("K0s") || sSpecies.EqualTo("Lambda")) {
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"gap00"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"CL1"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"PVz8"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"TPCcls90"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"V0sCPA099"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"V0sCrossFind1"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"V0sDaugDCA3"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"V0sDaugPt02"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"V0sDecRad1"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"V0sDecRad10"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            // FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"V0sFinderOn"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            // FlowFluct.Exec(Form("\"%s\",%d,\"%s\",%d,\"%s\",\"%s\",\"%s\"",sSpecies.Data(), iCent, Form("%s/%s/%s/",sPath.Data(),sSpecies.Data(),"V0sPVDCA3"), bCor, sHistoTwo.Data(), sHistoFour.Data(),sOutFileName.Data()));
            // } else if(sSpecies.EqualTo("Phi")) {
                //     single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"CL1",sBaseline.Data(),-1.0));
                //     single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"FB768",sBaseline.Data(),-1.0));
                //     single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PVz8",sBaseline.Data(),-1.0));
                //     single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"TPCcls90",sBaseline.Data(),-1.0));
                //     single.Exec(Form("\"%s\",\"%s\",\"%s\",\"%s\",\"%s\",%f",hist.Data(),path.Data(),pathOut.Data(),"PID3sigma",sBaseline.Data(),-1.0));
            } else {
                return kFALSE;
            }
    }

    return;
}

void FlowFluctAll()
{
    FlowFluctForSyst("Charged");
    FlowFluctForSyst("Pion");
    FlowFluctForSyst("Kaon");
    FlowFluctForSyst("Proton");
    FlowFluctForSyst("K0s");
    FlowFluctForSyst("Lambda");
}
