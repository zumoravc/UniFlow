AliAnalysisTaskChargedFlow* AddTaskChargedFlow(
								     UInt_t filterbit = 768,
								     Double_t etaCut = 0.8,
								     Double_t vtxCut = 10.,
								     Double_t minPt = 0.2,
								     Double_t maxPt = 5.0,
								     Short_t nHarm = 2,
                                     Bool_t IsSample = kTRUE,
								     Int_t noclus = 70,
                                     Int_t effsys = 0,
                                     Int_t charge = 0,
								     Bool_t IsLHC10h = kTRUE,
								     Bool_t IsPileUp = kTRUE,
								     Int_t nPtb = 14,
								     Short_t nCentFl = 0,
								     Bool_t ismc = kFALSE,
								     TString uniqueID = "",
										 Bool_t differentialFlow = kTRUE)
{
  // Creates a pid task and adds it to the analysis manager
  
  // Get the pointer to the existing analysis manager via the static
  //access method
  //=========================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskChargedFlow.C", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the
  // analysis manager The availability of MC handler can also be
  // checked here.
  // =========================================================================
  if (!mgr->GetInputEventHandler()) {
    Error("AddTaskChargedFlow", "This task requires an input event handler");
    return NULL;
  }  
  
  // Create the task and configure it 
  //========================================================================
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    
  Double_t ptBins[] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0};
  //Double_t ptBins[] = {0.2, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0};
    
  AliAnalysisTaskChargedFlow* taskFlowEp = new AliAnalysisTaskChargedFlow("taskFlowEp");
  taskFlowEp->SetDebugLevel(3);
  taskFlowEp->SetAnalysisType(type);
  taskFlowEp->SetFilterbit(filterbit);
  taskFlowEp->SetNoClus(noclus);
  taskFlowEp->SetCharge(charge);
  taskFlowEp->SetEff(effsys);
  taskFlowEp->SetEtaCut(etaCut);
  taskFlowEp->SetVtxCut(vtxCut);
  taskFlowEp->SetMinPt(minPt);
  taskFlowEp->SetMaxPt(maxPt);
  taskFlowEp->SetNHarmonic(nHarm);
  taskFlowEp->SetIsSample(IsSample);
  taskFlowEp->SetFlagLHC10h(IsLHC10h);
  taskFlowEp->SetFlagPileUp(IsPileUp);
  taskFlowEp->SetNPtBins(nPtb);
  taskFlowEp->SetPtBins(ptBins);
  taskFlowEp->SetAnalysisMC(ismc);
  taskFlowEp->SetCentFlag(nCentFl);
	taskFlowEp->SetDifferentialFlowFlag(differentialFlow);
  mgr->AddTask(taskFlowEp);
 
  
  // Create ONLY the output containers for the data produced by the
  // task.  Get and connect other common input/output containers via
  // the manager as below
   //=======================================================================
  //TString fileName = AliAnalysisManager::GetCommonFileName();
  //fileName+=suffixName;
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *cout_hist = mgr->CreateContainer(Form("output_%s", uniqueID.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, Form("v4New.root:%s", uniqueID.Data()));
  mgr->ConnectInput (taskFlowEp, 0, cinput);
  mgr->ConnectOutput(taskFlowEp, 1, cout_hist);
  
  // Return task pointer at the end
  return taskFlowEp;
}
