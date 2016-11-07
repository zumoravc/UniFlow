///////////////////////////////////////////////////////////////////
//                                                               //            
// AddFlowPID                                                     //
// Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2016       //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisTaskFlowPID* AddTaskFlowPID(TString name = "name")
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) return 0x0;
  if (!mgr->GetInputEventHandler())	return 0x0;

  TString fileName = AliAnalysisManager::GetCommonFileName();   // by default, a file is open for writing. here, we get the filename
  fileName += Form(":%s",name.Data());      // create a subfolder in the file
  
  AliAnalysisTaskFlowPID* task = new AliAnalysisTaskFlowPID(name.Data()); // now we create an instance of your task
  if(!task) return 0x0;

  mgr->AddTask(task); // add your task to the manager
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer()); // your task needs input: here we connect the manager to your task
  mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("Events_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data())); // same for the output
  mgr->ConnectOutput(task,2,mgr->CreateContainer(Form("Tracks_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data())); // same for the output
  mgr->ConnectOutput(task,3,mgr->CreateContainer(Form("PID_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data())); // same for the output
  mgr->ConnectOutput(task,4,mgr->CreateContainer(Form("V0s_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data())); // same for the output
  mgr->ConnectOutput(task,5,mgr->CreateContainer(Form("QA_%s",name.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data())); // same for the output
  return task;
}
