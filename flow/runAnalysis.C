void runAnalysis()
{
    Bool_t local = 0; // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t gridTest = 0; // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    
    TString sGridMode = "full";
    //TString sGridMode = "terminate";
    
    Bool_t bMergeViaJDL = kTRUE;
    //Bool_t bMergeViaJDL = kFALSE;

    TString sWorkDir = "V0s/13-QM12-check";
    TString sOutDir = "outFlow";
    
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
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskFlowPID");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    //Add this here: run before your task, but after definition of manager and input handler
    gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    AliMultSelectionTask* taskMultSelection = AddTaskMultSelection(kFALSE); // user mode:

    // PID response needed for V0s PID
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse* taskPIDResponse = AddTaskPIDResponse(kFALSE); // not MC

    gROOT->LoadMacro("AliAnalysisTaskFlowPID.cxx++g"); // compile the class (locally)
    gROOT->LoadMacro("AddTaskFlowPID.C"); // load the addtask macro
    
    AliAnalysisTaskFlowPID* taskFlowPID = AddTaskFlowPID("flowPID_JHEP"); // JHEP + V0 analysis notes
    // tracks & event selection cuts
    taskFlowPID->SetAODAnalysis(kTRUE);
    taskFlowPID->SetPbPbAnalysis(kTRUE);
    taskFlowPID->SetPeriod10h(kTRUE);
    taskFlowPID->SetCentFlag(0);
    taskFlowPID->SetPVtxZMax(10);
    taskFlowPID->SetTrackEtaMax(0.8);
    taskFlowPID->SetTrackPtMax(10.);
    taskFlowPID->SetTrackPtMin(0.1);
    taskFlowPID->SetNumTPCclsMin(70);
    taskFlowPID->SetTrackFilterBit(128);
    taskFlowPID->SetDiffFlow(kFALSE);
    taskFlowPID->SetPID(kTRUE);
    // V0 selection cuts
    taskFlowPID->SetV0sOnFly(kFALSE);
    taskFlowPID->SetV0sTPCRefit(kTRUE);
    taskFlowPID->SetV0sRejectKinks(kTRUE);
    taskFlowPID->SetV0sDCAPVMin(0.1);
    taskFlowPID->SetV0sDCAPVMax(0.);
    taskFlowPID->SetV0sDCADaughtersMax(1.);
    taskFlowPID->SetV0sDecayRadiusMin(5.);
    taskFlowPID->SetV0sDecayRadiusMax(0.);
    taskFlowPID->SetV0sDaughterPtMin(0.1);
    taskFlowPID->SetV0sDaughterEtaMax(0.8);
    taskFlowPID->SetV0sMotherEtaMax(0.8);
    taskFlowPID->SetV0sMotherRapMax(0.);
    taskFlowPID->SetV0sMotherPtMin(0.2);
    taskFlowPID->SetV0sMotherPtMax(5.);
    taskFlowPID->SetV0sK0sCPAMin(0.998);
    taskFlowPID->SetV0sLambdaCPAMin(0.998);
    taskFlowPID->SetV0sK0sNumTauMax(3.);
    taskFlowPID->SetV0sK0sArmenterosAlphaMin(0.2);
    taskFlowPID->SetV0sLambdaNumTauMax(3.);
    taskFlowPID->SetV0sProtonNumSigmaMax(3.);
    taskFlowPID->SetV0sProtonPIDPtMax(1.2);
    
    AliAnalysisTaskFlowPID* taskFlowPID_lose = AddTaskFlowPID("flowPID_lose"); // loser than JHEP
    // tracks & event selection cuts
    taskFlowPID_lose->SetAODAnalysis(kTRUE);
    taskFlowPID_lose->SetPbPbAnalysis(kTRUE);
    taskFlowPID_lose->SetPeriod10h(kTRUE);
    taskFlowPID_lose->SetCentFlag(0);
    taskFlowPID_lose->SetPVtxZMax(10.);
    taskFlowPID_lose->SetTrackEtaMax(0.8);
    taskFlowPID_lose->SetTrackPtMax(10.);
    taskFlowPID_lose->SetTrackPtMin(0.1);
    taskFlowPID_lose->SetNumTPCclsMin(70);
    taskFlowPID_lose->SetTrackFilterBit(128);
    taskFlowPID_lose->SetDiffFlow(0);
    taskFlowPID_lose->SetPID(kTRUE);
    // V0 selection cuts
    taskFlowPID_lose->SetV0sOnFly(kFALSE);
    taskFlowPID_lose->SetV0sTPCRefit(kTRUE);
    taskFlowPID_lose->SetV0sRejectKinks(kTRUE);
    taskFlowPID_lose->SetV0sDCAPVMin(0.05);
    taskFlowPID_lose->SetV0sDCAPVMax(0.);
    taskFlowPID_lose->SetV0sDCADaughtersMax(1.);
    taskFlowPID_lose->SetV0sDecayRadiusMin(3.);
    taskFlowPID_lose->SetV0sDecayRadiusMax(0.);
    taskFlowPID_lose->SetV0sDaughterPtMin(0.1);
    taskFlowPID_lose->SetV0sDaughterEtaMax(0.8);
    taskFlowPID_lose->SetV0sMotherEtaMax(0.8);
    taskFlowPID_lose->SetV0sMotherRapMax(0.);
    taskFlowPID_lose->SetV0sMotherPtMin(0.2);
    taskFlowPID_lose->SetV0sMotherPtMax(5.);
    taskFlowPID_lose->SetV0sK0sCPAMin(0.996);
    taskFlowPID_lose->SetV0sLambdaCPAMin(0.996);
    taskFlowPID_lose->SetV0sK0sArmenterosAlphaMin(0.2);
    taskFlowPID_lose->SetV0sK0sNumTauMax(4.);
    taskFlowPID_lose->SetV0sLambdaNumTauMax(4.);
    taskFlowPID_lose->SetV0sProtonNumSigmaMax(5.);
    taskFlowPID_lose->SetV0sProtonPIDPtMax(1.2);

    AliAnalysisTaskFlowPID* taskFlowPID_tight = AddTaskFlowPID("flowPID_tight"); // loser than JHEP
    // tracks & event selection cuts
    taskFlowPID_tight->SetAODAnalysis(kTRUE);
    taskFlowPID_tight->SetPbPbAnalysis(kTRUE);
    taskFlowPID_tight->SetPeriod10h(kTRUE);
    taskFlowPID_tight->SetCentFlag(0);
    taskFlowPID_tight->SetPVtxZMax(10);
    taskFlowPID_tight->SetTrackEtaMax(0.8);
    taskFlowPID_tight->SetTrackPtMax(10.);
    taskFlowPID_tight->SetTrackPtMin(0.1);
    taskFlowPID_tight->SetNumTPCclsMin(70);
    taskFlowPID_tight->SetTrackFilterBit(128);
    taskFlowPID_tight->SetDiffFlow(0);
    taskFlowPID_tight->SetPID(kTRUE);
    // V0 selection cuts
    taskFlowPID_tight->SetV0sOnFly(kFALSE);
    taskFlowPID_tight->SetV0sTPCRefit(kTRUE);
    taskFlowPID_tight->SetV0sRejectKinks(kTRUE);
    taskFlowPID_tight->SetV0sDCAPVMin(0.1);
    taskFlowPID_tight->SetV0sDCAPVMax(0.);
    taskFlowPID_tight->SetV0sDCADaughtersMax(1.);
    taskFlowPID_tight->SetV0sDecayRadiusMin(5.);
    taskFlowPID_tight->SetV0sDecayRadiusMax(0.);
    taskFlowPID_tight->SetV0sDaughterPtMin(0.1);
    taskFlowPID_tight->SetV0sDaughterEtaMax(0.8);
    taskFlowPID_tight->SetV0sMotherEtaMax(0.8);
    taskFlowPID_tight->SetV0sMotherRapMax(0.);
    taskFlowPID_tight->SetV0sMotherPtMin(0.2);
    taskFlowPID_tight->SetV0sMotherPtMax(5.);
    taskFlowPID_tight->SetV0sK0sCPAMin(0.998);
    taskFlowPID_tight->SetV0sLambdaCPAMin(0.999);
    taskFlowPID_tight->SetV0sK0sNumTauMax(2.5);
    taskFlowPID_tight->SetV0sK0sArmenterosAlphaMin(0.2);
    taskFlowPID_tight->SetV0sLambdaNumTauMax(2.5);
    taskFlowPID_tight->SetV0sProtonNumSigmaMax(2.5);
    taskFlowPID_tight->SetV0sProtonPIDPtMax(1.5);



    if (!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus(); 
    mgr->SetUseProgressBar(1, 100);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        chain->Add("~/NBI/Flow/flow/testData/2010/LHC10h/000138275/ESDs/pass2/AOD160/0803/AliAOD.root");  // add a few files to the chain (change this so that your local files are added)
        mgr->StartAnalysis("local", chain); // start the analysis locally, reading the events from the TChain
    } else {
        // if we want to run on grid, we create and configure the plugin
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        // also specify the include (header) paths on grid
        alienHandler->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        // make sure your source files get copied to grid
        alienHandler->SetAdditionalLibs("AliAnalysisTaskFlowPID.cxx AliAnalysisTaskFlowPID.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskFlowPID.cxx");
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
        
        // ++ 45 runs 
        // all
        //Int_t runNumber[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364};//..++
        // part1
        //Int_t runNumber[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871};
        //part2
        //Int_t runNumber[] = {138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364};//..++

        // -- 46 runs
        //Int_t runNumber[] = {138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};
        // part 1    
        //Int_t runNumber[] = {138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595};
        // part 2
        Int_t runNumber[] = {137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};

        Int_t iNumRuns = sizeof(runNumber) / sizeof(runNumber[0]);
        
        for (Int_t i = 0; i < iNumRuns; i++)
        {
            //if (i == sizeof(runArray) / sizeof(runArray[1])) break;
            alienHandler->AddRunNumber(runNumber[i]);
        }
        

        alienHandler->SetMasterResubmitThreshold(90);
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(200);
        alienHandler->SetExecutable("FlowPID.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(30000);
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