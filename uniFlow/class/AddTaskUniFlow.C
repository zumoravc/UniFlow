///////////////////////////////////////////////////////////////////
//
// AddTaskUniFlow.C macro
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016
//
//  See AliAnalysisTaskUniFlow(.cxx) for details & documentation
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;
class AliAnalysisTaskUniFlow;

AliAnalysisTaskUniFlow* AddTaskUniFlow(TString name, AliAnalysisTaskUniFlow::ColSystem colSys)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) { return NULL; }
  if (!mgr->GetInputEventHandler())	{ return NULL; }

  TString fileName = AliAnalysisManager::GetCommonFileName();   // by default, a file is open for writing. here, we get the filename
  fileName += Form(":%s",name.Data());      // create a subfolder in the file

  AliAnalysisTaskUniFlow* task = new AliAnalysisTaskUniFlow(name.Data(), colSys); // now we create an instance of your task
  if(!task) { return NULL; }

  // task default settings (ColSystem independent)
  task->SetAnalysisType(AliAnalysisTaskUniFlow::kAOD);
  task->SetRunMode(AliAnalysisTaskUniFlow::kFull);
  task->SetNumEventsAnalyse(50);
  task->SetTrigger(AliVEvent::kINT7);
  task->SetMC(kFALSE);
  task->SetFlowFillWeights(kTRUE);
  task->SetUseWeights3D(kFALSE);
  task->SetFillQAhistos(kTRUE);
  task->SetSampling(kFALSE);
  task->SetProcessPID(kTRUE);
  task->SetProcessPhi(kTRUE);
  task->SetProcessV0s(kTRUE);
  task->SetFlowRFPsPt(0.2,5.0);
  task->SetFlowPOIsPt(0.0,10.0);
  task->SetFlowEta(0.8);

  // task default settings dependent on ColSystem (colSys)
  if(colSys == AliAnalysisTaskUniFlow::ColSystem::kPbPb) {
    // Pb-Pb
    task->SetCentrality(AliAnalysisTaskUniFlow::kV0M);
    task->SetPVtxZMax(10.0);
    task->SetRejectAddPileUp(kFALSE);
    task->SetChargedNumTPCclsMin(70);
    task->SetChargedDCAzMax(0.0);
    task->SetChargedDCAxyMax(0.0);
    task->SetChargedTrackFilterBit(96);
    task->SetPIDUseAntiProtonOnly(kFALSE);
    task->SetPIDNumSigmasCombinedNoTOFrejection(kTRUE);
    task->SetPIDNumSigmasPionMax(0.0);
    task->SetPIDNumSigmasKaonMax(0.0);
    task->SetPIDNumSigmasProtonMax(0.0);
    task->SetPIDNumSigmasTPCRejectElectron(0.0);
    task->SetUseBayesPID(kTRUE);
    task->SetPIDBayesProbPionMin(0.95);
    task->SetPIDBayesProbKaonMin(0.85);
    task->SetPIDBayesProbProtonMin(0.85);
    task->SetV0sOnFly(kFALSE);
    task->SetV0sTPCRefit(kTRUE);
    task->SetV0sRejectKinks(kTRUE);
    task->SetV0sDaughterNumTPCClsMin(70);
    task->SetV0sDaughterNumTPCrossMin(70);
    task->SetV0sDaughterNumTPCFindMin(1);
    task->SetV0sDaughterNumTPCClsPIDMin(70);
    task->SetV0sDaughterRatioCrossFindMin(0.8);
    task->SetV0sUseCrossMassRejection(kTRUE);
    task->SetV0sCrossMassCutK0s(0.005);
    task->SetV0sCrossMassCutLambda(0.010);
    task->SetV0sDCAPVMin(0.1);
    task->SetV0sDCAPVMax(0.0);
    task->SetV0sDCAPVzMax(0.0);
    task->SetV0sDaughtersFilterBit(0);
    task->SetV0sDCADaughtersMin(0.0);
    task->SetV0sDCADaughtersMax(0.5);
    task->SetV0sDecayRadiusMin(5.0);
    task->SetV0sDecayRadiusMax(100.0);
    task->SetV0sDaughterEtaMax(0.8);
    task->SetV0sDaughterPtMin(0.0);
    task->SetV0sDaughterPtMax(0.0);
    task->SetV0sMotherRapMax(0.0);
    task->SetV0sK0sInvMassMin(0.4);
    task->SetV0sK0sInvMassMax(0.6);
    task->SetV0sLambdaInvMassMin(1.08);
    task->SetV0sLambdaInvMassMax(1.16);
    task->SetV0sK0sCPAMin(0.998);
    task->SetV0sLambdaCPAMin(0.998);
    task->SetV0sK0sNumTauMax(0.0);
    task->SetV0sLambdaNumTauMax(0.0);
    task->SetV0sK0sArmenterosAlphaMin(0.2);
    task->SetV0sLambdaArmenterosAlphaMax(0.0);
    task->SetV0sK0sPionNumTPCSigmaMax(3.0);
    task->SetV0sLambdaPionNumTPCSigmaMax(3.0);
    task->SetV0sLambdaProtonNumTPCSigmaMax(3.0);
    task->SetPhiInvMassMin(0.99);
    task->SetPhiInvMassMax(1.07);
  } else {
    // p-Pb & pp
    // NB: so far based on "previous" default ctor values
    task->SetCentrality(AliAnalysisTaskUniFlow::kV0A);
    task->SetPVtxZMax(10.0);
    task->SetRejectAddPileUp(kFALSE);
    task->SetChargedNumTPCclsMin(70);
    task->SetChargedDCAzMax(0.0);
    task->SetChargedDCAxyMax(0.0);
    task->SetChargedTrackFilterBit(96);
    task->SetPIDUseAntiProtonOnly(kFALSE);
    task->SetPIDNumSigmasCombinedNoTOFrejection(kTRUE);
    task->SetPIDNumSigmasPionMax(0.0);
    task->SetPIDNumSigmasKaonMax(0.0);
    task->SetPIDNumSigmasProtonMax(0.0);
    task->SetPIDNumSigmasTPCRejectElectron(0.0);
    task->SetUseBayesPID(kFALSE);
    task->SetPIDBayesProbPionMin(0.0);
    task->SetPIDBayesProbKaonMin(0.0);
    task->SetPIDBayesProbProtonMin(0.0);
    task->SetV0sOnFly(kFALSE);
    task->SetV0sTPCRefit(kTRUE);
    task->SetV0sRejectKinks(kTRUE);
    task->SetV0sDaughterNumTPCClsMin(0);
    task->SetV0sDaughterNumTPCrossMin(70);
    task->SetV0sDaughterNumTPCFindMin(1);
    task->SetV0sDaughterNumTPCClsPIDMin(0);
    task->SetV0sDaughterRatioCrossFindMin(0.8);
    task->SetV0sUseCrossMassRejection(kTRUE);
    task->SetV0sCrossMassCutK0s(0.005);
    task->SetV0sCrossMassCutLambda(0.010);
    task->SetV0sDCAPVMin(0.06);
    task->SetV0sDCAPVMax(0.0);
    task->SetV0sDCAPVzMax(0.0);
    task->SetV0sDaughtersFilterBit(0);
    task->SetV0sDCADaughtersMin(0.0);
    task->SetV0sDCADaughtersMax(1.0);
    task->SetV0sDecayRadiusMin(0.5);
    task->SetV0sDecayRadiusMax(200.0);
    task->SetV0sDaughterEtaMax(0.8);
    task->SetV0sDaughterPtMin(0.0);
    task->SetV0sDaughterPtMax(0.0);
    task->SetV0sMotherRapMax(0.0);
    task->SetV0sK0sInvMassMin(0.4);
    task->SetV0sK0sInvMassMax(0.6);
    task->SetV0sLambdaInvMassMin(1.08);
    task->SetV0sLambdaInvMassMax(1.16);
    task->SetV0sK0sCPAMin(0.97);
    task->SetV0sLambdaCPAMin(0.995);
    task->SetV0sK0sNumTauMax(7.46);
    task->SetV0sLambdaNumTauMax(3.8);
    task->SetV0sK0sArmenterosAlphaMin(0.2);
    task->SetV0sLambdaArmenterosAlphaMax(0.0);
    task->SetV0sK0sPionNumTPCSigmaMax(5.0);
    task->SetV0sLambdaPionNumTPCSigmaMax(5.0);
    task->SetV0sLambdaProtonNumTPCSigmaMax(5.0);
    task->SetPhiInvMassMin(0.99);
    task->SetPhiInvMassMax(1.07);
  }

  mgr->AddTask(task); // add your task to the manager

  // Creating containers
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer* cOutput1 = mgr->CreateContainer(Form("Flow_Refs_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput2 = mgr->CreateContainer(Form("Flow_Charged_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput3 = mgr->CreateContainer(Form("Flow_Pion_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput4 = mgr->CreateContainer(Form("Flow_Kaon_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput5 = mgr->CreateContainer(Form("Flow_Proton_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput6 = mgr->CreateContainer(Form("Flow_K0s_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput7 = mgr->CreateContainer(Form("Flow_Lambda_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput8 = mgr->CreateContainer(Form("Flow_Phi_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput9 = mgr->CreateContainer(Form("QA_Events_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput10 = mgr->CreateContainer(Form("QA_Charged_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput11 = mgr->CreateContainer(Form("QA_PID_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput12 = mgr->CreateContainer(Form("QA_V0s_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput13 = mgr->CreateContainer(Form("QA_Phi_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput14 = mgr->CreateContainer(Form("Flow_Weights_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));

  // Connecting containers to task
  mgr->ConnectInput(task,0,cInput0); // your task needs input: here we connect the manager to your task
  mgr->ConnectOutput(task,1,cOutput1);
  mgr->ConnectOutput(task,2,cOutput2);
  mgr->ConnectOutput(task,3,cOutput3);
  mgr->ConnectOutput(task,4,cOutput4);
  mgr->ConnectOutput(task,5,cOutput5);
  mgr->ConnectOutput(task,6,cOutput6);
  mgr->ConnectOutput(task,7,cOutput7);
  mgr->ConnectOutput(task,8,cOutput8);
  mgr->ConnectOutput(task,9,cOutput9);
  mgr->ConnectOutput(task,10,cOutput10);
  mgr->ConnectOutput(task,11,cOutput11);
  mgr->ConnectOutput(task,12,cOutput12);
  mgr->ConnectOutput(task,13,cOutput13);
  mgr->ConnectOutput(task,14,cOutput14);

  return task;
}
