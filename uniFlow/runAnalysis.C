void runAnalysis()
{
    Bool_t local = 1; // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t gridTest = 0; // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)

    TString sGridMode = "full";
    //TString sGridMode = "terminate";

    Bool_t bMergeViaJDL = kTRUE;
    //Bool_t bMergeViaJDL = kFALSE;

    TString sWorkDir = "";
    TString sOutDir = "outFlow";

    // run switcher
    // ++ 45 runs
        // all
        //Int_t runNumber[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364};//..++
        // part1
        Int_t runNumber[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871};
        //part2
        //Int_t runNumber[] = {138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364};//..++
        // testing sample
        //Int_t runNumber[] = {138870, 138837, 138732, 138730, 138666, 138662};


    // -- 46 runs
        //Int_t runNumber[] = {138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};
        // part 1
        //Int_t runNumber[] = {138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595};
        // part 2
        //Int_t runNumber[] = {137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};




    // since we will compile a class, tell root where to look for headers
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");

    gSystem->AddIncludePath("-I$ALICE_ROOT/include");
    gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
    //gSystem->AddIncludePath("-I$ALICE_PHYSICS/../src/OADB/COMMON/MULTIPLICITY");
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
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskUniFlow");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    // Physics selection as suggested by HMTF
    gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE,kTRUE);

    //Add this here: run before your task, but after definition of manager and input handler
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask* taskMultSelection = AddTaskMultSelection(kFALSE); // user mode:
    taskMultSelection->SetSelectedTriggerClass(AliVEvent::kMB);

    // PID response needed for PID
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse* taskPIDResponse = AddTaskPIDResponse(kFALSE); // not MC

    /*
    // PID response QA by ALICE
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa* taskPIDqa = AddTaskPIDqa("PID_QA.root");
    */

    gROOT->LoadMacro("AliAnalysisTaskUniFlow.cxx++g"); // compile the class (locally)
    gROOT->LoadMacro("AddTaskUniFlow.C"); // load the addtask macro

    AliAnalysisTaskUniFlow* task1 = AddTaskUniFlow("UniFlow_test");
    // Analysis
    task1->SetRunMode(AliAnalysisTaskUniFlow::kFull);
    task1->SetNumEventsAnalyse(10);
    task1->SetAnalysisType(AliAnalysisTaskUniFlow::kAOD);
    task1->SetSampling(kFALSE);
    task1->SetProcessCharged(kTRUE);
    // task1->SetProcessPID(kTRUE);
    task1->SetProcessV0s(kTRUE);
    // Flow
    task1->SetFlowRFPsPtMin(0.2);
    task1->SetFlowRFPsPtMax(5.);
    task1->SetFlowChargedPtMin(0.2);
    task1->SetFlowChargedPtMax(5.);
    // Events selection
    task1->SetTrigger(0);
    task1->SetColisionSystem(AliAnalysisTaskUniFlow::kPP);
    task1->SetPeriod(AliAnalysisTaskUniFlow::k16l);
    task1->SetPVtxZMax(10);
    // Charged selection
    task1->SetChargedEtaMax(0.8);
    // task1->SetChargedPtMin(0.2);
    // task1->SetChargedPtMax(5.);
    // task1->SetChargedDCAzMax(0.1);
    // task1->SetChargedDCAxyMax(0.2);
    task1->SetChargedNumTPCclsMin(70);
    task1->SetChargedTrackFilterBit(96);
    // PID selection
    task1->SetPionNumSigmasMax(3);
    task1->SetKaonNumSigmasMax(3);
    task1->SetProtonNumSigmasMax(3);
    task1->SetPIDBayesProbPionMin(0.95);
    task1->SetPIDBayesProbKaonMin(0.8);
    task1->SetPIDBayesProbProtonMin(0.8);
    // V0 selection cuts
    task1->SetV0sOnFly(kFALSE);
    task1->SetV0sTPCRefit(kTRUE);
    task1->SetV0sRejectKinks(kTRUE);
    task1->SetV0sUseCrossMassRejection(kTRUE);
    task1->SetV0sDCAPVMin(0.1);
    task1->SetV0sDCAPVMax(0.);
    task1->SetV0sDCADaughtersMax(1.);
    task1->SetV0sDecayRadiusMin(5.);
    task1->SetV0sDecayRadiusMax(0.);
    task1->SetV0sDaughterPtMin(0.1);
    task1->SetV0sDaughterEtaMax(0.8);
    task1->SetV0sMotherEtaMax(0.8);
    task1->SetV0sMotherRapMax(0.);
    task1->SetV0sMotherPtMin(0.2);
    task1->SetV0sMotherPtMax(10.);
    task1->SetV0sK0sCPAMin(0.998);
    task1->SetV0sLambdaCPAMin(0.998);
    task1->SetV0sK0sNumTauMax(3.);
    task1->SetV0sK0sArmenterosAlphaMin(0.2);
    task1->SetV0sLambdaNumTauMax(3.);
    task1->SetV0sProtonNumSigmaMax(3.);
    task1->SetV0sProtonPIDPtMax(1.2);

    /*
    AliAnalysisTaskUniFlow* task2 = AddTaskUniFlow("UniFlow_FB768_Nsigma");
    // tracks & event selection cuts
    task2->SetAODAnalysis(kTRUE);
    task2->SetPbPbAnalysis(kTRUE);
    task2->SetPPbAnalysis(kFALSE);
    task2->SetPPAnalysis(kFALSE);
    task2->SetPeriod10h(kFALSE);
    task2->SetCentFlag(0);
    task2->SetPVtxZMax(10);
    task2->SetTrackEtaMax(0.8);
    task2->SetTrackPtMin(0.2);
    task2->SetTrackPtMax(5.);
    task2->SetNumTPCclsMin(70);
    task2->SetTrackFilterBit(768);
    task2->SetPionNumSigmasMax(3);
    task2->SetKaonNumSigmasMax(3);
    task2->SetProtonNumSigmasMax(3);
    //task2->SetDoFlow(kTRUE);
    task2->SetDiffFlow(kTRUE);
    task2->SetPID(kTRUE);
    task2->SetDoV0s(kTRUE);
    task2->SetSampling(kTRUE);
    task2->SetDoFlowGenFramKatarina(kTRUE);
    task2->SetDoOldFlow(kFALSE);
    task2->SetUseBayesPID(kFALSE);
    task2->SetPIDBayesProbPionMin(0.95);
    task2->SetPIDBayesProbKaonMin(0.95);
    task2->SetPIDBayesProbProtonMin(0.95);
    // V0 selection cuts
    task2->SetV0sOnFly(kFALSE);
    task2->SetV0sTPCRefit(kTRUE);
    task2->SetV0sRejectKinks(kTRUE);
    task2->SetV0sDCAPVMin(0.1);
    task2->SetV0sDCAPVMax(0.);
    task2->SetV0sDCADaughtersMax(1.);
    task2->SetV0sDecayRadiusMin(5.);
    task2->SetV0sDecayRadiusMax(0.);
    task2->SetV0sDaughterPtMin(0.1);
    task2->SetV0sDaughterEtaMax(0.8);
    task2->SetV0sMotherEtaMax(0.8);
    task2->SetV0sMotherRapMax(0.);
    task2->SetV0sMotherPtMin(0.2);
    task2->SetV0sMotherPtMax(10.);
    task2->SetV0sK0sCPAMin(0.998);
    task2->SetV0sLambdaCPAMin(0.998);
    task2->SetV0sK0sNumTauMax(3.);
    task2->SetV0sK0sArmenterosAlphaMin(0.2);
    task2->SetV0sLambdaNumTauMax(3.);
    task2->SetV0sProtonNumSigmaMax(3.);
    task2->SetV0sProtonPIDPtMax(1.2);
    */
    if (!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);
    //mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 500);

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
        alienHandler->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskUniFlow.cxx AliAnalysisTaskUniFlow.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskUniFlow.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20160907-1");
        //alienHandler->SetAliPhysicsVersion("vAN-20160131-1");
        // select the input data
        alienHandler->SetGridDataDir("/alice/data/2010/LHC10h/");
        alienHandler->SetDataPattern("ESDs/pass2/AOD160/*/AliAOD.root");
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
        alienHandler->SetSplitMaxInputFileNumber(200);
        alienHandler->SetExecutable("UniFlow.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(36000);
        alienHandler->SetJDLName("UniFlow.jdl");
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
