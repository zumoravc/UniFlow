#!/bin/bash

# ===========================================================
# Script for generating runPP.C to given folder
# Arguments:
#   1 - output
#		2 - tag
#   3 - period
# ===========================================================

echo "##### Generating runPPbRun2.sh #####"

# checking parameters
if [ $# -ne 3 ]; then
	echo "Wrong number of parameters: $# (3 needed). Exit!"
	exit
fi

if [ "$1" = "" ] || [ "$2" = "" ] ||  "$3" = "" ]; then
	echo "Empty parameter given. Exit!"
	exit
fi

output=${1}
tag=${2}
period="$(tr '[:lower:]' '[:upper:]' <<< ${3:0:3})${3:3}"

if [ ! -d "${output}" ]; then
	#path folder does not exist
  echo "Output folder (arg. 1) does not exists. Exit!"
	exit
fi

# if [[ ! "${period}" == "LHC16k"* ]] && [ ! "${period}" == "LHC16l"* ]]; then
# 	# runlist file does not exist
#   echo "Wrong period '${period}' (not LHC16k nor LHC16l). Exit!"
# 	exit
# fi
# arguments passed
echo "--- Listing arguments ---"
echo "output:\"${output}\""
echo "tag:\"${tag}\""
echo "period:\"${period}\""
echo "-------------------------"

# writing following script into $file
cat > ${output}/runPPbRun2.C <<EOT
void runPPbRun2()
{
    Bool_t local = 1; // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t gridTest = 0; // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)

    TString sGridMode = "full";
    //TString sGridMode = "terminate";

    Bool_t bMergeViaJDL = kTRUE;
    //Bool_t bMergeViaJDL = kFALSE;

		TString sWorkDir = "${tag}";
    TString sOutDir = "output";
    TString sPeriod = "${period}";

    // run switcher
		// Run2 5.02 TeV
		// RunList_LHC16r_pass1_CentralBarrelTracking_hadronPID_20170202_v0.txt [12 runs]
    Int_t runNumber[] = {266318, 266317, 266316, 266208, 266197, 266196, 266187, 265754, 265744, 265607, 265596, 265594};
    // RunList_LHC16s_pass1_CentralBarrelTracking_hadronPID_20170202_v0.txt [16 runs]
    // Int_t runNumber[] = {267110, 267081, 267077, 267072, 267070, 267030, 266998, 266997, 266994, 266993, 266944, 266886, 266885, 266883, 266882, 266437};

		// Run2 8.16 TeV
		// RunList_LHC16t_pass1_CentralBarrelTracking_hadronPID_20170202_v0.txt [4 runs]
		// Int_t runNumber[] = {267166, 267165, 267164, 267163};
		// RunList_LHC16q_pass1_CentralBarrelTracking_hadronPID_20170202_v0.txt [32 runs]
		// Int_t runNumber[] = {265525, 265521, 265501, 265500, 265499, 265435, 265427, 265426, 265425, 265424, 265422, 265421, 265420, 265419, 265388, 265387, 265385, 265384, 265383, 265381, 265378, 265377, 265344, 265343, 265342, 265339, 265338, 265336, 265335, 265334, 265332, 265309};

    // since we will compile a class, tell root where to look for headers
    gROOT->ProcessLine(".include \$ROOTSYS/include");
    gROOT->ProcessLine(".include \$ALICE_ROOT/include");

    gSystem->AddIncludePath("-I\$ALICE_ROOT/include");
    gSystem->AddIncludePath("-I\$ALICE_PHYSICS/include");
    //gSystem->AddIncludePath("-I\$ALICE_PHYSICS/../src/OADB/COMMON/MULTIPLICITY");
    gSystem->Load("libPhysics.so");
    gSystem->Load("libOADB.so");
    gSystem->Load("libSTEERBase");
    gSystem->Load("libESD");
    gSystem->Load("libAOD");
    gSystem->Load("libANALYSIS");
    gSystem->Load("libANALYSISalice");
    gSystem->Load("libANALYSISaliceBase");
    gSystem->Load("libCORRFW");

    // create the analysis manager
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskFlowPID");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    //Add this here: run before your task, but after definition of manager and input handler

    //gROOT->LoadMacro("\$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    //AliMultSelectionTask* taskMultSelection = AddTaskMultSelection(kFALSE); // user mode:
    //taskMultSelection->SetSelectedTriggerClass(AliVEvent::kMB);

    // Physics selection as suggested by HMTF
    gROOT->ProcessLine(".L \$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE,kTRUE);

    // PID response needed for PID
    gROOT->LoadMacro("\$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse* taskPIDResponse = AddTaskPIDResponse(kFALSE); // not MC

    /*
    // PID response QA by ALICE
    gROOT->LoadMacro("\$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa* taskPIDqa = AddTaskPIDqa("PID_QA.root");
    */

    gROOT->LoadMacro("AliAnalysisTaskFlowPID.cxx++g"); // compile the class (locally)
    gROOT->LoadMacro("AddTaskFlowPID.C"); // load the addtask macro


    AliAnalysisTaskFlowPID* task1 = AddTaskFlowPID("flowPID_CENT_wSDD");
    //task1->SelectCollisionCandidates(AliVEvent::kINT7)
    // analysis cuts & switches
    task1->SetAODAnalysis(kTRUE);
    task1->SetPPAnalysis(kTRUE);
    task1->SetTrigger(0); // kINT7
    task1->SetSampling(kTRUE);
    //task1->SetUseOldCent(kFALSE);
    //task1->SetCentFlag(0);
    //task1->SetPeriod10h(kFALSE);
    //task1->SetRejectPileUpSPD(kFALSE);
    task1->SetDoFlow(kTRUE);
    task1->SetDoFlowGenFramKatarina(kTRUE);
    task1->SetDoOldFlow(kFALSE);
    task1->SetDiffFlow(kTRUE);
    task1->SetTracksScan(kFALSE);
    task1->SetPID(kTRUE);
    task1->SetDoV0s(kTRUE);
    // event selection
    task1->SetPVtxZMax(10);
    task1->SetTrackEtaMax(0.8);
    task1->SetTrackPtMin(0.2);
    task1->SetTrackPtMax(5.);
    task1->SetTrackDCAzMax(0.0);
    task1->SetNumTPCclsMin(70);
    task1->SetTrackFilterBit(96);
    task1->SetUseBayesPID(kTRUE);
    task1->SetPIDBayesProbPionMin(0.95);
    task1->SetPIDBayesProbKaonMin(0.8);
    task1->SetPIDBayesProbProtonMin(0.8);
    //task1->SetPionNumSigmasMax(3);
    //task1->SetKaonNumSigmasMax(3);
    //task1->SetProtonNumSigmasMax(3);
    // V0 selection cut
    task1->SetV0sOnFly(kFALSE);
    task1->SetV0sTPCRefit(kTRUE);
    task1->SetV0sRejectKinks(kTRUE);
    task1->SetV0sDCAPVMin(0.06);
    task1->SetV0sDCAPVMax(0.);
    task1->SetV0sDCADaughtersMax(1.);
    task1->SetV0sDecayRadiusMin(0.5);
    task1->SetV0sDecayRadiusMax(0.);
    task1->SetV0sDaughterPtMin(0.1);
    task1->SetV0sDaughterEtaMax(0.8);
    task1->SetV0sMotherEtaMax(0.8);
    task1->SetV0sMotherRapMax(0.);
    task1->SetV0sMotherPtMin(0.2);
    task1->SetV0sMotherPtMax(5.);
    task1->SetV0sK0sCPAMin(0.97);
    task1->SetV0sLambdaCPAMin(0.995);
    task1->SetV0sK0sNumTauMax(5.);
    task1->SetV0sK0sArmenterosAlphaMin(0.2);
    task1->SetV0sLambdaNumTauMax(3.8);
    task1->SetV0sProtonNumSigmaMax(0);
    task1->SetV0sProtonPIDPtMax(0);

    // AliAnalysisTaskFlowPID* task2 = AddTaskFlowPID("flowPID_kINT7_noSampled");
    // //task2->SelectCollisionCandidates(AliVEvent::kINT7)
    // // analysis cuts & switches
    // task2->SetAODAnalysis(kTRUE);
    // task2->SetPPAnalysis(kTRUE);
    // task2->SetTrigger(0); // kINT7
    // task2->SetSampling(kFALSE);
    // //task2->SetUseOldCent(kFALSE);
    // //task2->SetCentFlag(0);
    // //task2->SetPeriod10h(kFALSE);
    // //task2->SetRejectPileUpSPD(kFALSE);
    // task2->SetDoFlow(kTRUE);
    // task2->SetDoFlowGenFramKatarina(kTRUE);
    // task2->SetDoOldFlow(kFALSE);
    // task2->SetDiffFlow(kTRUE);
    // task2->SetTracksScan(kFALSE);
    // task2->SetPID(kFALSE);
    // task2->SetDoV0s(kFALSE);
    // // event selection
    // task2->SetPVtxZMax(10);
    // task2->SetTrackEtaMax(0.8);
    // task2->SetTrackPtMin(0.2);
    // task2->SetTrackPtMax(5.);
    // task2->SetTrackDCAzMax(0.0);
    // task2->SetNumTPCclsMin(70);
    // task2->SetTrackFilterBit(96);
    // task2->SetUseBayesPID(kTRUE);
    // task2->SetPIDBayesProbPionMin(0.8);
    // task2->SetPIDBayesProbKaonMin(0.8);
    // task2->SetPIDBayesProbProtonMin(0.8);
    // //task2->SetPionNumSigmasMax(3);
    // //task2->SetKaonNumSigmasMax(3);
    // //task2->SetProtonNumSigmasMax(3);
    // // V0 selection cut
    // task2->SetV0sOnFly(kFALSE);
    // task2->SetV0sTPCRefit(kTRUE);
    // task2->SetV0sRejectKinks(kTRUE);
    // task2->SetV0sDCAPVMin(0.1);
    // task2->SetV0sDCAPVMax(0.);
    // task2->SetV0sDCADaughtersMax(1.);
    // task2->SetV0sDecayRadiusMin(5.);
    // task2->SetV0sDecayRadiusMax(0.);
    // task2->SetV0sDaughterPtMin(0.1);
    // task2->SetV0sDaughterEtaMax(0.8);
    // task2->SetV0sMotherEtaMax(0.8);
    // task2->SetV0sMotherRapMax(0.);
    // task2->SetV0sMotherPtMin(0.2);
    // task2->SetV0sMotherPtMax(10.);
    // task2->SetV0sK0sCPAMin(0.998);
    // task2->SetV0sLambdaCPAMin(0.998);
    // task2->SetV0sK0sNumTauMax(3.);
    // task2->SetV0sK0sArmenterosAlphaMin(0.2);
    // task2->SetV0sLambdaNumTauMax(3.);
    // task2->SetV0sProtonNumSigmaMax(3.);
    // task2->SetV0sProtonPIDPtMax(1.2);


    if (!mgr->InitAnalysis()) return;
    //git mgr->SetDebugLevel(2);
    //mgr->PrintStatus();
    //mgr->SetUseProgressBar(1, 100);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        //chain->Add("~/NBI/Flow/flow/testData/2010/LHC10h/000138275/ESDs/pass2/AOD160/0803/AliAOD.root");  // add a few files to the chain (change this so that your local files are added)
        chain->Add("~/NBI/Flow/flow/testData/2016/LHC16l/000259888/pass1/AOD/015/AliAOD.root");  // add a few files to the chain (change this so that your local files are added)
        mgr->StartAnalysis("local", chain); // start the analysis locally, reading the events from the TChain
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I\$ROOTSYS/include -I\$ALICE_ROOT/include -I\$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskFlowPID.cxx AliAnalysisTaskFlowPID.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskFlowPID.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20161219-1");
        // select the input data
        alienHandler->SetGridDataDir(Form("/alice/data/2016/%s/",sPeriod.Data()));
        alienHandler->SetDataPattern("/pass1_CENT_wSDD/AOD/*/AliAOD.root");

        // MC has no prefix, data has prefix 000
        alienHandler->SetRunPrefix("000");
        // runnumber

        Int_t iNumRuns = sizeof(runNumber) / sizeof(runNumber[0]);

        for (Int_t i = 0; i < iNumRuns; i++)
        {
            //if (i == sizeof(runArray) / sizeof(runArray[1])) break;
            alienHandler->AddRunNumber(runNumber[i]);
        }


        alienHandler->SetMasterResubmitThreshold(90);
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(15);
        alienHandler->SetExecutable("FlowPID.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(50000);
        alienHandler->SetJDLName("FlowPID.jdl");
        alienHandler->SetPrice(1);
        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);

        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(bMergeViaJDL);

        // define the output folders
        alienHandler->SetGridWorkingDir(sWorkDir.Data());
        alienHandler->SetGridOutputDir(sOutDir.Data());

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode(sGridMode.Data());
            mgr->StartAnalysis("grid");
        }
    }
}
EOT

if [ -f ${output}/runPPbRun2.C ]; then
	echo "File generated!"
else
	echo "File NOT generated!"
	exit
fi
