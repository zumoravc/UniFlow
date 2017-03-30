///////////////////////////////////////////////////////////////////
//                                                               //
// AddUniFlow                                                     //
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016       //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskUniFlow* AddTaskUniFlow(TString name = "name")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;

  TString fileName = AliAnalysisManager::GetCommonFileName();   // by default, a file is open for writing. here, we get the filename
  fileName += Form(":%s",name.Data());      // create a subfolder in the file

  AliAnalysisTaskUniFlow* task = new AliAnalysisTaskUniFlow(name.Data()); // now we create an instance of your task
  if(!task) return 0x0;

  mgr->AddTask(task); // add your task to the manager

  // Creating containers
  AliAnalysisDataContainer* cInput0 = mgr->GetCommonInputContainer();
  // AliAnalysisDataContainer* cOutput1 = mgr->CreateContainer(Form("Cumulants_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput1 = mgr->CreateContainer(Form("Events_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput2 = mgr->CreateContainer(Form("Charged_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput3 = mgr->CreateContainer(Form("PID_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  AliAnalysisDataContainer* cOutput4 = mgr->CreateContainer(Form("V0s_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  // AliAnalysisDataContainer* cOutput4 = mgr->CreateContainer(Form("PID_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  // AliAnalysisDataContainer* cOutput5 = mgr->CreateContainer(Form("V0s_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));
  //AliAnalysisDataContainer* cOutput6 = mgr->CreateContainer(Form("QA_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("%s",fileName.Data()));

  // Connecting containers to task
  mgr->ConnectInput(task,0,cInput0); // your task needs input: here we connect the manager to your task
  mgr->ConnectOutput(task,1,cOutput1);
  mgr->ConnectOutput(task,2,cOutput2);
  mgr->ConnectOutput(task,3,cOutput3);
  mgr->ConnectOutput(task,4,cOutput4);
  // mgr->ConnectOutput(task,2,cOutput2);
  // mgr->ConnectOutput(task,3,cOutput3);
  // mgr->ConnectOutput(task,4,cOutput4);
  // mgr->ConnectOutput(task,5,cOutput5);
  //mgr->ConnectOutput(task,6,cOutput6);

  return task;
}
