void runPP()
{
    Bool_t local = 0; // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t gridTest = 0; // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    
    TString sGridMode = "full";
    //TString sGridMode = "terminate";
    
    Bool_t bMergeViaJDL = kTRUE;
    //Bool_t bMergeViaJDL = kFALSE;

    TString sWorkDir = "pp-2016k-QA2";
    TString sOutDir = "outFlow";
    //TString sPeriod = "LHC16l";
    TString sPeriod = "LHC16k";

    // run switcher
    // RunList_LHC16l_pass1_CentralBarrelTracking_hadronPID_20161122_v1.txt [75 runs]
    //Int_t runNumber[] = {260014, 260011, 260010, 259888, 259868, 259867, 259866, 259860, 259842, 259841, 259822, 259789, 259788, 259781, 259756, 259752, 259751, 259750, 259748, 259747, 259713, 259711, 259705, 259704, 259703, 259700, 259697, 259668, 259650, 259649, 259477, 259473, 259396, 259395, 259394, 259389, 259388, 259382, 259378, 259342, 259341, 259340, 259339, 259336, 259334, 259307, 259305, 259303, 259302, 259274, 259273, 259272, 259271, 259270, 259269, 259264, 259263, 259261, 259257, 259204, 259164, 259162, 259118, 259117, 259099, 259096, 259091, 259090, 259088, 258964, 258962, 258923, 258921, 258920, 258919};
    
    // RunList_LHC16k_pass1_CentralBarrelTracking_hadronPID_20161121_v0.txt [97 runs]
    // part 1 [44 runs]
    //Int_t runNumber[] = {258537, 258499, 258477, 258456, 258454, 258452, 258426, 258393, 258391, 258387, 258359, 258336, 258332, 258307, 258306, 258303, 258302, 258301, 258299, 258278, 258274, 258273, 258271, 258270, 258258, 258257, 258256, 258204, 258203, 258202, 258198, 258197, 258178, 258117, 258114, 258113, 258109, 258108, 258107, 258063, 258062, 258060, 258059, 258053, 258049, 258045, 258042, 258041, 258039};
    // part 2 [43 runs]
    Int_t runNumber[] = {258019, 258017, 258014, 258012, 258008, 258003, 257992, 257989, 257986, /*257979,*/ 257963, 257960, 257957, 257939, 257937, 257936, 257892, 257855, 257853, 257851, 257850, 257804, 257803, 257800, 257799, 257798, 257797, 257773, 257765, 257757, 257754, 257737, 257735, 257734, 257733, 257727, 257725, 257724, 257697, 257694, 257692, 257691, 257689, 257688, 257687, 257685, 257684, 257682};
    
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

    //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
    //AliMultSelectionTask* taskMultSelection = AddTaskMultSelection(kFALSE); // user mode:
    //taskMultSelection->SetSelectedTriggerClass(AliVEvent::kMB);
    
    // Physics selection as suggested by HMTF
    gROOT->ProcessLine(".L $ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
    AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(kFALSE,kTRUE);

    // PID response needed for PID
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
    AliAnalysisTaskPIDResponse* taskPIDResponse = AddTaskPIDResponse(kFALSE); // not MC

    /*
    // PID response QA by ALICE
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDqa.C");
    AliAnalysisTaskPIDqa* taskPIDqa = AddTaskPIDqa("PID_QA.root");
    */

    gROOT->LoadMacro("AliAnalysisTaskFlowPID.cxx++g"); // compile the class (locally)
    gROOT->LoadMacro("AddTaskFlowPID.C"); // load the addtask macro
    


    AliAnalysisTaskFlowPID* task1 = AddTaskFlowPID("flowPID_FB768_kINT7");
    //task1->SelectCollisionCandidates(AliVEvent::kINT7)
    // analysis cuts & switches
    task1->SetAODAnalysis(kTRUE);
    task1->SetPPAnalysis(kTRUE);
    task1->SetTrigger(0); // kINT7
    task1->SetSampling(kFALSE);
    //task1->SetUseOldCent(kFALSE);
    //task1->SetCentFlag(0);
    //task1->SetPeriod10h(kFALSE);
    //task1->SetRejectPileUpSPD(kFALSE);
    task1->SetDoFlowGenFramKatarina(kTRUE);
    task1->SetDoOldFlow(kFALSE);
    //task1->SetDoFlow(kTRUE);
    task1->SetDiffFlow(kTRUE);
    task1->SetPID(kTRUE);
    task1->SetDoV0s(kTRUE);
    // event selection
    task1->SetPVtxZMax(10);
    task1->SetTrackEtaMax(0.8);
    task1->SetTrackPtMin(0.2);
    task1->SetTrackPtMax(5.);
    task1->SetNumTPCclsMin(70);
    task1->SetTrackFilterBit(768);
    task1->SetUseBayesPID(kTRUE);
    task1->SetPIDBayesProbPionMin(0.8);
    task1->SetPIDBayesProbKaonMin(0.8);
    task1->SetPIDBayesProbProtonMin(0.8);
    //task1->SetPionNumSigmasMax(3);
    //task1->SetKaonNumSigmasMax(3);
    //task1->SetProtonNumSigmasMax(3);
    // V0 selection cuts
    task1->SetV0sOnFly(kFALSE);
    task1->SetV0sTPCRefit(kTRUE);
    task1->SetV0sRejectKinks(kTRUE);
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
        
    AliAnalysisTaskFlowPID* task2 = AddTaskFlowPID("flowPID_FB768_kHighMultV0");
    //task2->SelectCollisionCandidates(AliVEvent::kINT7)
    // analysis cuts & switches
    task2->SetAODAnalysis(kTRUE);
    task2->SetPPAnalysis(kTRUE);
    task2->SetTrigger(1); // kINT7
    task2->SetSampling(kFALSE);
    //task2->SetUseOldCent(kFALSE);
    //task2->SetCentFlag(0);
    //task2->SetPeriod10h(kFALSE);
    //task2->SetRejectPileUpSPD(kFALSE);
    task2->SetDoFlowGenFramKatarina(kTRUE);
    task2->SetDoOldFlow(kFALSE);
    //task2->SetDoFlow(kTRUE);
    task2->SetDiffFlow(kTRUE);
    task2->SetPID(kTRUE);
    task2->SetDoV0s(kTRUE);
    // event selection
    task2->SetPVtxZMax(10);
    task2->SetTrackEtaMax(0.8);
    task2->SetTrackPtMin(0.2);
    task2->SetTrackPtMax(5.);
    task2->SetNumTPCclsMin(70);
    task2->SetTrackFilterBit(768);
    task2->SetUseBayesPID(kTRUE);
    task2->SetPIDBayesProbPionMin(0.8);
    task2->SetPIDBayesProbKaonMin(0.8);
    task2->SetPIDBayesProbProtonMin(0.8);
    //task2->SetPionNumSigmasMax(3);
    //task2->SetKaonNumSigmasMax(3);
    //task2->SetProtonNumSigmasMax(3);
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
    
    AliAnalysisTaskFlowPID* task3 = AddTaskFlowPID("flowPID_FB1_kINT7");
    //task3->SelectCollisionCandidates(AliVEvent::kINT7)
    // analysis cuts & switches
    task3->SetAODAnalysis(kTRUE);
    task3->SetPPAnalysis(kTRUE);
    task3->SetTrigger(0); // kINT7
    task3->SetSampling(kFALSE);
    //task3->SetUseOldCent(kFALSE);
    //task3->SetCentFlag(0);
    //task3->SetPeriod10h(kFALSE);
    //task3->SetRejectPileUpSPD(kFALSE);
    task3->SetDoFlowGenFramKatarina(kTRUE);
    task3->SetDoOldFlow(kFALSE);
    //task3->SetDoFlow(kTRUE);
    task3->SetDiffFlow(kTRUE);
    task3->SetPID(kTRUE);
    task3->SetDoV0s(kTRUE);
    // event selection
    task3->SetPVtxZMax(10);
    task3->SetTrackEtaMax(0.8);
    task3->SetTrackPtMin(0.2);
    task3->SetTrackPtMax(5.);
    task3->SetNumTPCclsMin(70);
    task3->SetTrackFilterBit(1);
    task3->SetUseBayesPID(kTRUE);
    task3->SetPIDBayesProbPionMin(0.8);
    task3->SetPIDBayesProbKaonMin(0.8);
    task3->SetPIDBayesProbProtonMin(0.8);
    //task3->SetPionNumSigmasMax(3);
    //task3->SetKaonNumSigmasMax(3);
    //task3->SetProtonNumSigmasMax(3);
    // V0 selection cuts
    task3->SetV0sOnFly(kFALSE);
    task3->SetV0sTPCRefit(kTRUE);
    task3->SetV0sRejectKinks(kTRUE);
    task3->SetV0sDCAPVMin(0.1);
    task3->SetV0sDCAPVMax(0.);
    task3->SetV0sDCADaughtersMax(1.);
    task3->SetV0sDecayRadiusMin(5.);
    task3->SetV0sDecayRadiusMax(0.);
    task3->SetV0sDaughterPtMin(0.1);
    task3->SetV0sDaughterEtaMax(0.8);
    task3->SetV0sMotherEtaMax(0.8);
    task3->SetV0sMotherRapMax(0.);
    task3->SetV0sMotherPtMin(0.2);
    task3->SetV0sMotherPtMax(10.);
    task3->SetV0sK0sCPAMin(0.998);
    task3->SetV0sLambdaCPAMin(0.998);
    task3->SetV0sK0sNumTauMax(3.);
    task3->SetV0sK0sArmenterosAlphaMin(0.2);
    task3->SetV0sLambdaNumTauMax(3.);
    task3->SetV0sProtonNumSigmaMax(3.);
    task3->SetV0sProtonPIDPtMax(1.2);  

    AliAnalysisTaskFlowPID* task4 = AddTaskFlowPID("flowPID_FB1_kHighMultV0");
    //task4->SelectCollisionCandidates(AliVEvent::kINT7)
    // analysis cuts & switches
    task4->SetAODAnalysis(kTRUE);
    task4->SetPPAnalysis(kTRUE);
    task4->SetTrigger(1); // kINT7
    task4->SetSampling(kFALSE);
    //task4->SetUseOldCent(kFALSE);
    //task4->SetCentFlag(0);
    //task4->SetPeriod10h(kFALSE);
    //task4->SetRejectPileUpSPD(kFALSE);
    task4->SetDoFlowGenFramKatarina(kTRUE);
    task4->SetDoOldFlow(kFALSE);
    //task4->SetDoFlow(kTRUE);
    task4->SetDiffFlow(kTRUE);
    task4->SetPID(kTRUE);
    task4->SetDoV0s(kTRUE);
    // event selection
    task4->SetPVtxZMax(10);
    task4->SetTrackEtaMax(0.8);
    task4->SetTrackPtMin(0.2);
    task4->SetTrackPtMax(5.);
    task4->SetNumTPCclsMin(70);
    task4->SetTrackFilterBit(1);
    task4->SetUseBayesPID(kTRUE);
    task4->SetPIDBayesProbPionMin(0.8);
    task4->SetPIDBayesProbKaonMin(0.8);
    task4->SetPIDBayesProbProtonMin(0.8);
    //task4->SetPionNumSigmasMax(3);
    //task4->SetKaonNumSigmasMax(3);
    //task4->SetProtonNumSigmasMax(3);
    // V0 selection cuts
    task4->SetV0sOnFly(kFALSE);
    task4->SetV0sTPCRefit(kTRUE);
    task4->SetV0sRejectKinks(kTRUE);
    task4->SetV0sDCAPVMin(0.1);
    task4->SetV0sDCAPVMax(0.);
    task4->SetV0sDCADaughtersMax(1.);
    task4->SetV0sDecayRadiusMin(5.);
    task4->SetV0sDecayRadiusMax(0.);
    task4->SetV0sDaughterPtMin(0.1);
    task4->SetV0sDaughterEtaMax(0.8);
    task4->SetV0sMotherEtaMax(0.8);
    task4->SetV0sMotherRapMax(0.);
    task4->SetV0sMotherPtMin(0.2);
    task4->SetV0sMotherPtMax(10.);
    task4->SetV0sK0sCPAMin(0.998);
    task4->SetV0sLambdaCPAMin(0.998);
    task4->SetV0sK0sNumTauMax(3.);
    task4->SetV0sK0sArmenterosAlphaMin(0.2);
    task4->SetV0sLambdaNumTauMax(3.);
    task4->SetV0sProtonNumSigmaMax(3.);
    task4->SetV0sProtonPIDPtMax(1.2); 
    /*
    AliAnalysisTaskFlowPID* task5 = AddTaskFlowPID("flowPID_FB128_kINT7");
    //task5->SelectCollisionCandidates(AliVEvent::kINT7)
    // analysis cuts & switches
    task5->SetAODAnalysis(kTRUE);
    task5->SetPPAnalysis(kTRUE);
    task5->SetTrigger(0); // kINT7
    task5->SetSampling(kFALSE);
    //task5->SetUseOldCent(kFALSE);
    //task5->SetCentFlag(0);
    //task5->SetPeriod10h(kFALSE);
    //task5->SetRejectPileUpSPD(kFALSE);
    task5->SetDoFlowGenFramKatarina(kTRUE);
    task5->SetDoOldFlow(kFALSE);
    //task5->SetDoFlow(kTRUE);
    task5->SetDiffFlow(kTRUE);
    task5->SetPID(kTRUE);
    task5->SetDoV0s(kTRUE);
    // event selection
    task5->SetPVtxZMax(10);
    task5->SetTrackEtaMax(0.8);
    task5->SetTrackPtMin(0.2);
    task5->SetTrackPtMax(5.);
    task5->SetNumTPCclsMin(70);
    task5->SetTrackFilterBit(128);
    task5->SetUseBayesPID(kTRUE);
    task5->SetPIDBayesProbPionMin(0.8);
    task5->SetPIDBayesProbKaonMin(0.8);
    task5->SetPIDBayesProbProtonMin(0.8);
    //task5->SetPionNumSigmasMax(3);
    //task5->SetKaonNumSigmasMax(3);
    //task5->SetProtonNumSigmasMax(3);
    // V0 selection cuts
    task5->SetV0sOnFly(kFALSE);
    task5->SetV0sTPCRefit(kTRUE);
    task5->SetV0sRejectKinks(kTRUE);
    task5->SetV0sDCAPVMin(0.1);
    task5->SetV0sDCAPVMax(0.);
    task5->SetV0sDCADaughtersMax(1.);
    task5->SetV0sDecayRadiusMin(5.);
    task5->SetV0sDecayRadiusMax(0.);
    task5->SetV0sDaughterPtMin(0.1);
    task5->SetV0sDaughterEtaMax(0.8);
    task5->SetV0sMotherEtaMax(0.8);
    task5->SetV0sMotherRapMax(0.);
    task5->SetV0sMotherPtMin(0.2);
    task5->SetV0sMotherPtMax(10.);
    task5->SetV0sK0sCPAMin(0.998);
    task5->SetV0sLambdaCPAMin(0.998);
    task5->SetV0sK0sNumTauMax(3.);
    task5->SetV0sK0sArmenterosAlphaMin(0.2);
    task5->SetV0sLambdaNumTauMax(3.);
    task5->SetV0sProtonNumSigmaMax(3.);
    task5->SetV0sProtonPIDPtMax(1.2);  

    AliAnalysisTaskFlowPID* task6 = AddTaskFlowPID("flowPID_FB128_kHighMultV0");
    //task6->SelectCollisionCandidates(AliVEvent::kINT7)
    // analysis cuts & switches
    task6->SetAODAnalysis(kTRUE);
    task6->SetPPAnalysis(kTRUE);
    task6->SetTrigger(1); // kINT7
    task6->SetSampling(kFALSE);
    //task6->SetUseOldCent(kFALSE);
    //task6->SetCentFlag(0);
    //task6->SetPeriod10h(kFALSE);
    //task6->SetRejectPileUpSPD(kFALSE);
    task6->SetDoFlowGenFramKatarina(kTRUE);
    task6->SetDoOldFlow(kFALSE);
    //task6->SetDoFlow(kTRUE);
    task6->SetDiffFlow(kTRUE);
    task6->SetPID(kTRUE);
    task6->SetDoV0s(kTRUE);
    // event selection
    task6->SetPVtxZMax(10);
    task6->SetTrackEtaMax(0.8);
    task6->SetTrackPtMin(0.2);
    task6->SetTrackPtMax(5.);
    task6->SetNumTPCclsMin(70);
    task6->SetTrackFilterBit(128);
    task6->SetUseBayesPID(kTRUE);
    task6->SetPIDBayesProbPionMin(0.8);
    task6->SetPIDBayesProbKaonMin(0.8);
    task6->SetPIDBayesProbProtonMin(0.8);
    //task6->SetPionNumSigmasMax(3);
    //task6->SetKaonNumSigmasMax(3);
    //task6->SetProtonNumSigmasMax(3);
    // V0 selection cuts
    task6->SetV0sOnFly(kFALSE);
    task6->SetV0sTPCRefit(kTRUE);
    task6->SetV0sRejectKinks(kTRUE);
    task6->SetV0sDCAPVMin(0.1);
    task6->SetV0sDCAPVMax(0.);
    task6->SetV0sDCADaughtersMax(1.);
    task6->SetV0sDecayRadiusMin(5.);
    task6->SetV0sDecayRadiusMax(0.);
    task6->SetV0sDaughterPtMin(0.1);
    task6->SetV0sDaughterEtaMax(0.8);
    task6->SetV0sMotherEtaMax(0.8);
    task6->SetV0sMotherRapMax(0.);
    task6->SetV0sMotherPtMin(0.2);
    task6->SetV0sMotherPtMax(10.);
    task6->SetV0sK0sCPAMin(0.998);
    task6->SetV0sLambdaCPAMin(0.998);
    task6->SetV0sK0sNumTauMax(3.);
    task6->SetV0sK0sArmenterosAlphaMin(0.2);
    task6->SetV0sLambdaNumTauMax(3.);
    task6->SetV0sProtonNumSigmaMax(3.);
    task6->SetV0sProtonPIDPtMax(1.2); 
    */
    
    if (!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus(); 
    mgr->SetUseProgressBar(1, 100);

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
        alienHandler->SetAdditionalLibs("AliAnalysisTaskFlowPID.cxx AliAnalysisTaskFlowPID.h");
        alienHandler->SetAnalysisSource("AliAnalysisTaskFlowPID.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20161130-1");
        //alienHandler->SetAliPhysicsVersion("vAN-20160907-1");
        //alienHandler->SetAliPhysicsVersion("vAN-20160131-1");
        // select the input data
        alienHandler->SetGridDataDir(Form("/alice/data/2016/%s/",sPeriod.Data()));
        alienHandler->SetDataPattern("/pass1/AOD/*/AliAOD.root"); 

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
        alienHandler->SetSplitMaxInputFileNumber(100);
        alienHandler->SetExecutable("FlowPID.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(36000);
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