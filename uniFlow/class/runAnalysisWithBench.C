void runAnalysisWithBench()
{
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
    taskMultSelection->SetSelectedTriggerClass(AliVEvent::kINT7);

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

    AliAnalysisTaskUniFlow* task1 = AddTaskUniFlow("UniFlow");
    // Analysis
    task1->SetRunMode(AliAnalysisTaskUniFlow::kFull);
    // task1->SetRunMode(AliAnalysisTaskUniFlow::kTest);
    task1->SetNumEventsAnalyse(50);
    task1->SetAnalysisType(AliAnalysisTaskUniFlow::kAOD);
    // task1->SetSampling(kTRUE);
    // task1->SetFillQAhistos(kFALSE);
    task1->SetProcessCharged(kTRUE);
    task1->SetProcessPID(kTRUE);
    task1->SetProcessPhi(kTRUE);
    task1->SetProcessV0s(kTRUE);
    // Flow
    task1->SetFlowRFPsPtMin(0.2);
    task1->SetFlowRFPsPtMax(5.);
    // task1->SetFlowDoFourCorrelations(kFALSE);
    task1->SetFlowFillWeights(kFALSE);
    // task1->SetUseWeigthsFile("alien:///alice/cern.ch/user/v/vpacik/weights_preliminary_16q.root",kFALSE);
    // Events selection
    // task1->SetUseAliEventCuts();
    task1->SetTrigger(0);
    task1->SetColisionSystem(AliAnalysisTaskUniFlow::kPPb);
    task1->SetMultEstimator("V0A");
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
    task1->SetPIDUseAntiProtonOnly(kFALSE);
    task1->SetPIDNumSigmasPionMax(3);
    task1->SetPIDNumSigmasKaonMax(3);
    task1->SetPIDNumSigmasProtonMax(3);
    task1->SetUseBayesPID(kTRUE);
    task1->SetPIDBayesProbPionMin(0.95);
    task1->SetPIDBayesProbKaonMin(0.85);
    task1->SetPIDBayesProbProtonMin(0.85);
    // Phi selection
    // task1->SetPhiMotherEtaMax(0.8);
    // V0 selection cuts
    task1->SetV0sOnFly(kFALSE);
    task1->SetV0sTPCRefit(kTRUE);
    task1->SetV0sRejectKinks(kTRUE);
    task1->SetV0sUseCrossMassRejection(kTRUE);
    task1->SetV0sCrossMassCutK0s(0.005);
    task1->SetV0sCrossMassCutLambda(0.020);
    task1->SetV0sDCAPVMin(0.06);
    task1->SetV0sDCAPVMax(0.);
    // task1->SetV0sDCAPVzMax(1.);
    // task1->SetV0sDaughtersFilterBit(211);
    task1->SetV0sDCADaughtersMax(1.);
    task1->SetV0sDecayRadiusMin(0.5);
    task1->SetV0sDecayRadiusMax(0.);
    task1->SetV0sDaughterPtMin(0.1);
    task1->SetV0sDaughterEtaMax(0.8);
    task1->SetV0sMotherEtaMax(0.8);
    task1->SetV0sMotherRapMax(0.);
    task1->SetV0sK0sCPAMin(0.97);
    task1->SetV0sLambdaCPAMin(0.99);
    task1->SetV0sK0sNumTauMax(5);
    task1->SetV0sK0sArmenterosAlphaMin(0.2);
    task1->SetV0sLambdaNumTauMax(3.8);
    // task1->SetV0sK0sKaonNumTPCSigmaMax(3.);
    // task1->SetV0sLambdaPionNumTPCSigmaMax(3.);
    // task1->SetV0sLambdaProtonNumTPCSigmaMax(3.);


    if (!mgr->InitAnalysis()) return;
    //mgr->SetDebugLevel(2);
    //mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 500);

    // if you want to run locally, we need to define some input
    TChain* chain = new TChain("aodTree");
    //chain->Add("~/NBI/Flow/flow/testData/2010/LHC10h/000138275/ESDs/pass2/AOD160/0803/AliAOD.root");
    // chain->Add("~/NBI/Flow/flow/testData/2016/LHC16l/000259888/pass1/AOD/015/AliAOD.root");
    chain->Add("~/NBI/Flow/data/2016/LHC16q/000265427/pass1_CENT_wSDD/AOD/001/AliAOD.root");

    Int_t iNumRuns = 20;

    TH1D* hCPU = new TH1D("hCPU","CPU time dist; CPU time; Counts", 1000, 0, 100);
    TH1D* hReal = new TH1D("hReal","Real time dist; Real time; Counts", 1000, 0, 100);

    Double_t dCPUtime = 0.0;
    Double_t dRealtime = 0.0;

    TStopwatch* watch = new TStopwatch();
    for(Int_t iRun(0); iRun < iNumRuns; ++iRun)
    {
      watch->Reset();

      printf("Starting analysis\n");
      watch->Start();
      mgr->StartAnalysis("local", chain); // start the analysis locally, reading the events from the TChain
      watch->Stop();
      dCPUtime = watch->CpuTime();
      dRealtime = watch->RealTime();
      printf("Stopped analysis\n");
      printf("%d :\t CPU %f\t Real %f\n",iRun,dCPUtime, dRealtime);
      hCPU->Fill(dCPUtime);
      hReal->Fill(dRealtime);
    }

    TCanvas* canBench = new TCanvas("canBench","canBench",1600,800);
    canBench->Divide(2,1);
    canBench->cd(1);
    hCPU->Draw();
    canBench->cd(2);
    hReal->Draw();

    canBench->SaveAs("bench_time_new.pdf");
}
