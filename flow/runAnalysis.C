void runAnalysis()
{
    Bool_t local = 1; // set if you want to run the analysis locally (kTRUE), or on grid (kFALSE)
    Bool_t gridTest = 0; // if you run on grid, specify test mode (kTRUE) or full grid model (kFALSE)
    
    TString sGridMode = "full";
    //TString sGridMode = "terminate";
    
    Bool_t bMergeViaJDL = kTRUE;
    //Bool_t bMergeViaJDL = kFALSE;

    TString sWorkDir = "11-GFK-PbPb-R1";
    TString sOutDir = "outFlow";
    
    // run switcher
    // ++ 45 runs 
        // all
        //Int_t runNumber[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364};//..++
        // part1
        //Int_t runNumber[] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871};
        //part2
        //Int_t runNumber[] = {138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364};//..++
        // testing sample
        //Int_t runNumber[] = {138870, 138837, 138732, 138730, 138666, 138662};


    // -- 46 runs
        //Int_t runNumber[] = {138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};
        // part 1    
        Int_t runNumber[] = {138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595};
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
    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskFlowPID");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

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

    gROOT->LoadMacro("AliAnalysisTaskFlowPID.cxx++g"); // compile the class (locally)
    gROOT->LoadMacro("AddTaskFlowPID.C"); // load the addtask macro
    
    AliAnalysisTaskFlowPID* task1 = AddTaskFlowPID("flowPID_FB768_New_PileOFF_PeriodOFF");
    // tracks & event selection cuts
    task1->SetAODAnalysis(kTRUE);
    task1->SetPbPbAnalysis(kFALSE);
    task1->SetPPAnalysis(kFALSE);
    task1->SetPPbAnalysis(kFALSE);
    task1->SetUseOldCent(kFALSE);
    task1->SetPeriod10h(kFALSE);
    task1->SetRejectPileUpSPD(kFALSE);
    task1->SetDoOldFlow(kFALSE);
    task1->SetDoFlowGenFramKatarina(kTRUE);
    task1->SetDiffFlow(kTRUE);
    task1->SetPID(kFALSE);
    task1->SetDoV0s(kFALSE);
    task1->SetSampling(kFALSE);
    task1->SetCentFlag(0);
    task1->SetPVtxZMax(10);
    task1->SetTrackEtaMax(0.8);
    task1->SetTrackPtMin(0.2);
    task1->SetTrackPtMax(5.);
    task1->SetNumTPCclsMin(70);
    task1->SetTrackFilterBit(768);
    task1->SetPionNumSigmasMax(3);
    task1->SetKaonNumSigmasMax(3);
    task1->SetProtonNumSigmasMax(3);
    //task1->SetDoFlow(kTRUE);
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

    /*
    AliAnalysisTaskFlowPID* task2 = AddTaskFlowPID("flowPID_FB768_NewCent_NoSampling");
    // tracks & event selection cuts
    task2->SetAODAnalysis(kTRUE);
    task2->SetPbPbAnalysis(kTRUE);
    task2->SetUseOldCent(kFALSE);
    task2->SetPeriod10h(kFALSE);
    task2->SetRejectPileUpSPD(kFALSE);
    task2->SetDoOldFlow(kFALSE);
    task2->SetDoFlowGenFramKatarina(kTRUE);
    task2->SetDiffFlow(kTRUE);
    task2->SetPID(kFALSE);
    task2->SetDoV0s(kFALSE);
    task2->SetSampling(kFALSE);
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
    
    AliAnalysisTaskFlowPID* task3 = AddTaskFlowPID("flowPID_FB768_OldCent_NoSampling");
    // tracks & event selection cuts
    task3->SetAODAnalysis(kTRUE);
    task3->SetPbPbAnalysis(kTRUE);
    task3->SetUseOldCent(kTRUE);
    task3->SetPeriod10h(kFALSE);
    task3->SetRejectPileUpSPD(kFALSE);
    task3->SetDoOldFlow(kFALSE);
    task3->SetDoFlowGenFramKatarina(kTRUE);
    task3->SetDiffFlow(kTRUE);
    task3->SetPID(kFALSE);
    task3->SetDoV0s(kFALSE);
    task3->SetSampling(kFALSE);
    task3->SetCentFlag(0);
    task3->SetPVtxZMax(10);
    task3->SetTrackEtaMax(0.8);
    task3->SetTrackPtMin(0.2);
    task3->SetTrackPtMax(5.);
    task3->SetNumTPCclsMin(70);
    task3->SetTrackFilterBit(768);
    task3->SetPionNumSigmasMax(3);
    task3->SetKaonNumSigmasMax(3);
    task3->SetProtonNumSigmasMax(3);
    //task3->SetDoFlow(kTRUE);
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
    
    */

    /*
     // Testing PileRejection & period ON/OFF
    AliAnalysisTaskFlowPID* task2 = AddTaskFlowPID("flowPID_FB768_New_PileON_PeriodOFF");
    // tracks & event selection cuts
    task2->SetAODAnalysis(kTRUE);
    task2->SetPbPbAnalysis(kTRUE);
    task2->SetPeriod10h(kFALSE);
    task2->SetRejectPileUpSPD(kTRUE);
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
    task2->SetPID(kFALSE);
    task2->SetDoV0s(kTRUE);
    task2->SetSampling(kFALSE);
    task2->SetDoFlowGenFramKatarina(kTRUE);
    task2->SetDoOldFlow(kFALSE);
    task2->SetUseOldCent(kFALSE);

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

    AliAnalysisTaskFlowPID* task3 = AddTaskFlowPID("flowPID_FB768_New_PileOFF_PeriodON");
    // tracks & event selection cuts
    task3->SetAODAnalysis(kTRUE);
    task3->SetPbPbAnalysis(kTRUE);
    task3->SetPeriod10h(kTRUE);
    task3->SetRejectPileUpSPD(kFALSE);
    task3->SetCentFlag(0);
    task3->SetPVtxZMax(10);
    task3->SetTrackEtaMax(0.8);
    task3->SetTrackPtMin(0.2);
    task3->SetTrackPtMax(5.);
    task3->SetNumTPCclsMin(70);
    task3->SetTrackFilterBit(768);
    task3->SetPionNumSigmasMax(3);
    task3->SetKaonNumSigmasMax(3);
    task3->SetProtonNumSigmasMax(3);
    //task3->SetDoFlow(kTRUE);
    task3->SetDiffFlow(kTRUE);
    task3->SetPID(kFALSE);
    task3->SetDoV0s(kTRUE);
    task3->SetSampling(kFALSE);
    task3->SetDoFlowGenFramKatarina(kTRUE);
    task3->SetDoOldFlow(kFALSE);
    task3->SetUseOldCent(kFALSE);

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

    AliAnalysisTaskFlowPID* task4 = AddTaskFlowPID("flowPID_FB768_New_PileON_PeriodON");
    // tracks & event selection cuts
    task4->SetAODAnalysis(kTRUE);
    task4->SetPbPbAnalysis(kTRUE);
    task4->SetPeriod10h(kTRUE);
    task4->SetRejectPileUpSPD(kTRUE);
    task4->SetCentFlag(0);
    task4->SetPVtxZMax(10);
    task4->SetTrackEtaMax(0.8);
    task4->SetTrackPtMin(0.2);
    task4->SetTrackPtMax(5.);
    task4->SetNumTPCclsMin(70);
    task4->SetTrackFilterBit(768);
    task4->SetPionNumSigmasMax(3);
    task4->SetKaonNumSigmasMax(3);
    task4->SetProtonNumSigmasMax(3);
    //task4->SetDoFlow(kTRUE);
    task4->SetDiffFlow(kTRUE);
    task4->SetPID(kFALSE);
    task4->SetDoV0s(kTRUE);
    task4->SetSampling(kFALSE);
    task4->SetDoFlowGenFramKatarina(kTRUE);
    task4->SetDoOldFlow(kFALSE);
    task4->SetUseOldCent(kFALSE);

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

    AliAnalysisTaskFlowPID* task8 = AddTaskFlowPID("flowPID_FB768_Old_PileON_PeriodON");
    // tracks & event selection cuts
    task8->SetAODAnalysis(kTRUE);
    task8->SetPbPbAnalysis(kTRUE);
    task8->SetPeriod10h(kTRUE);
    task8->SetRejectPileUpSPD(kTRUE);
    task8->SetCentFlag(0);
    task8->SetPVtxZMax(10);
    task8->SetTrackEtaMax(0.8);
    task8->SetTrackPtMin(0.2);
    task8->SetTrackPtMax(5.);
    task8->SetNumTPCclsMin(70);
    task8->SetTrackFilterBit(768);
    task8->SetPionNumSigmasMax(3);
    task8->SetKaonNumSigmasMax(3);
    task8->SetProtonNumSigmasMax(3);
    //task8->SetDoFlow(kTRUE);
    task8->SetDiffFlow(kTRUE);
    task8->SetPID(kFALSE);
    task8->SetDoV0s(kTRUE);
    task8->SetSampling(kFALSE);
    task8->SetDoFlowGenFramKatarina(kTRUE);
    task8->SetDoOldFlow(kFALSE);
    task8->SetUseOldCent(kTRUE);

    // V0 selection cuts
    task8->SetV0sOnFly(kFALSE);
    task8->SetV0sTPCRefit(kTRUE);
    task8->SetV0sRejectKinks(kTRUE);
    task8->SetV0sDCAPVMin(0.1);
    task8->SetV0sDCAPVMax(0.);
    task8->SetV0sDCADaughtersMax(1.);
    task8->SetV0sDecayRadiusMin(5.);
    task8->SetV0sDecayRadiusMax(0.);
    task8->SetV0sDaughterPtMin(0.1);
    task8->SetV0sDaughterEtaMax(0.8);
    task8->SetV0sMotherEtaMax(0.8);
    task8->SetV0sMotherRapMax(0.);
    task8->SetV0sMotherPtMin(0.2);
    task8->SetV0sMotherPtMax(10.);
    task8->SetV0sK0sCPAMin(0.998);
    task8->SetV0sLambdaCPAMin(0.998);
    task8->SetV0sK0sNumTauMax(3.);
    task8->SetV0sK0sArmenterosAlphaMin(0.2);
    task8->SetV0sLambdaNumTauMax(3.);
    task8->SetV0sProtonNumSigmaMax(3.);
    task8->SetV0sProtonPIDPtMax(1.2);

    AliAnalysisTaskFlowPID* task5 = AddTaskFlowPID("flowPID_FB768_Old_PileON_PeriodOFF");
    // tracks & event selection cuts
    task5->SetAODAnalysis(kTRUE);
    task5->SetPbPbAnalysis(kTRUE);
    task5->SetPeriod10h(kFALSE);
    task5->SetRejectPileUpSPD(kTRUE);
    task5->SetCentFlag(0);
    task5->SetPVtxZMax(10);
    task5->SetTrackEtaMax(0.8);
    task5->SetTrackPtMin(0.2);
    task5->SetTrackPtMax(5.);
    task5->SetNumTPCclsMin(70);
    task5->SetTrackFilterBit(768);
    task5->SetPionNumSigmasMax(3);
    task5->SetKaonNumSigmasMax(3);
    task5->SetProtonNumSigmasMax(3);
    //task5->SetDoFlow(kTRUE);
    task5->SetDiffFlow(kTRUE);
    task5->SetPID(kFALSE);
    task5->SetDoV0s(kTRUE);
    task5->SetSampling(kFALSE);
    task5->SetDoFlowGenFramKatarina(kTRUE);
    task5->SetDoOldFlow(kFALSE);
    task5->SetUseOldCent(kTRUE);

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

    AliAnalysisTaskFlowPID* task6 = AddTaskFlowPID("flowPID_FB768_Old_PileOFF_PeriodON");
    // tracks & event selection cuts
    task6->SetAODAnalysis(kTRUE);
    task6->SetPbPbAnalysis(kTRUE);
    task6->SetPeriod10h(kTRUE);
    task6->SetRejectPileUpSPD(kFALSE);
    task6->SetCentFlag(0);
    task6->SetPVtxZMax(10);
    task6->SetTrackEtaMax(0.8);
    task6->SetTrackPtMin(0.2);
    task6->SetTrackPtMax(5.);
    task6->SetNumTPCclsMin(70);
    task6->SetTrackFilterBit(768);
    task6->SetPionNumSigmasMax(3);
    task6->SetKaonNumSigmasMax(3);
    task6->SetProtonNumSigmasMax(3);
    //task6->SetDoFlow(kTRUE);
    task6->SetDiffFlow(kTRUE);
    task6->SetPID(kFALSE);
    task6->SetDoV0s(kTRUE);
    task6->SetSampling(kFALSE);
    task6->SetDoFlowGenFramKatarina(kTRUE);
    task6->SetDoOldFlow(kFALSE);
    task6->SetUseOldCent(kTRUE);

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

    AliAnalysisTaskFlowPID* task7 = AddTaskFlowPID("flowPID_FB768_Old_PileOFF_PeriodOFF");
    // tracks & event selection cuts
    task7->SetAODAnalysis(kTRUE);
    task7->SetPbPbAnalysis(kTRUE);
    task7->SetPeriod10h(kFALSE);
    task7->SetRejectPileUpSPD(kFALSE);
    task7->SetCentFlag(0);
    task7->SetPVtxZMax(10);
    task7->SetTrackEtaMax(0.8);
    task7->SetTrackPtMin(0.2);
    task7->SetTrackPtMax(5.);
    task7->SetNumTPCclsMin(70);
    task7->SetTrackFilterBit(768);
    task7->SetPionNumSigmasMax(3);
    task7->SetKaonNumSigmasMax(3);
    task7->SetProtonNumSigmasMax(3);
    //task7->SetDoFlow(kTRUE);
    task7->SetDiffFlow(kTRUE);
    task7->SetPID(kFALSE);
    task7->SetDoV0s(kTRUE);
    task7->SetSampling(kFALSE);
    task7->SetDoFlowGenFramKatarina(kTRUE);
    task7->SetDoOldFlow(kFALSE);
    task7->SetUseOldCent(kTRUE);

    // V0 selection cuts
    task7->SetV0sOnFly(kFALSE);
    task7->SetV0sTPCRefit(kTRUE);
    task7->SetV0sRejectKinks(kTRUE);
    task7->SetV0sDCAPVMin(0.1);
    task7->SetV0sDCAPVMax(0.);
    task7->SetV0sDCADaughtersMax(1.);
    task7->SetV0sDecayRadiusMin(5.);
    task7->SetV0sDecayRadiusMax(0.);
    task7->SetV0sDaughterPtMin(0.1);
    task7->SetV0sDaughterEtaMax(0.8);
    task7->SetV0sMotherEtaMax(0.8);
    task7->SetV0sMotherRapMax(0.);
    task7->SetV0sMotherPtMin(0.2);
    task7->SetV0sMotherPtMax(10.);
    task7->SetV0sK0sCPAMin(0.998);
    task7->SetV0sLambdaCPAMin(0.998);
    task7->SetV0sK0sNumTauMax(3.);
    task7->SetV0sK0sArmenterosAlphaMin(0.2);
    task7->SetV0sLambdaNumTauMax(3.);
    task7->SetV0sProtonNumSigmaMax(3.);
    task7->SetV0sProtonPIDPtMax(1.2);
    
    
    */

    AliAnalysisTaskFlowPID* task2 = AddTaskFlowPID("flowPID_FB768_Nsigma");
    // tracks & event selection cuts
    task2->SetAODAnalysis(kTRUE);
    task2->SetPbPbAnalysis(kTRUE);
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

    if (!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus(); 
    mgr->SetUseProgressBar(1, 100);

    if(local) {
        // if you want to run locally, we need to define some input
        TChain* chain = new TChain("aodTree");
        //chain->Add("~/NBI/Flow/flow/testData/2010/LHC10h/000138275/ESDs/pass2/AOD160/0803/AliAOD.root");  // add a few files to the chain (change this so that your local files are added)
        chain->Add("~/NBI/Flow/flow/testData/2010/LHC10h/000138870/ESDs/pass2/AOD160/0058/AliAOD.root");  // add a few files to the chain (change this so that your local files are added)
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