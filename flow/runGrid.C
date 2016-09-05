void runGrid()
{
  // loading libraries
  // since we will compile a class, tell root where to look for headers  
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");

  // Load common libraries
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");

  // Use AliRoot includes to compile our task
  //  gROOT->ProcessLine(".include $ALICE_ROOT/include"); 
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/../src/OADB/COMMON/MULTIPLICITY"); 
  gSystem->Load("libPhysics.so");
  gSystem->Load("libOADB.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libANALYSISaliceBase");
  gSystem->Load("libCORRFW");


  // create analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("flowEP");
  AliAODInputHandler *aodH = new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  // create an instance of your analysis task
  gROOT->LoadMacro("AliAnalysisTaskChargedFlow.cxx+g");
  gROOT->LoadMacro("AddTaskChargedFlow.C");


  // ============== running on grid ===================
  // if we want to run on grid, we create and configure the alienHandler
  AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
  // also specify the include (header) paths on grid
  alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
  // make sure your source files get copied to grid
  alienHandler->SetAdditionalLibs("AliAnalysisTaskChargedFlow.cxx AliAnalysisTaskChargedFlow.h");
  alienHandler->SetAnalysisSource("AliAnalysisTaskChargedFlow.cxx");
  // select the aliphysics version. All other packages are LOADED AUTOMATICALLY!
  alienHandler->SetAliPhysicsVersion("vAN-20160131-1");
  //alienHandler->SetAPIVersion("V1.1x"); // ???
  // select the input data
  alienHandler->SetGridDataDir("/alice/data/2015/LHC15o/");
  alienHandler->SetDataPattern("/pass_lowint_firstphys/AOD/*/AliAOD.root"); //..for ++ polarity);
  // MC has no prefix, data has prefix 000
  alienHandler->SetRunPrefix("000");
  // runnumber
  alienHandler->AddRunNumber(244918);
  // number of files per subjob
  alienHandler->SetSplitMaxInputFileNumber(40);
//  alienHandler->SetExecutable("myTask.sh");
  // specify how many seconds your job may take
  alienHandler->SetTTL(10000);
//  alienHandler->SetJDLName("myTask.jdl");

  alienHandler->SetOutputToRunNo(kTRUE);
  alienHandler->SetKeepLogs(kTRUE);
  // merging: run with kTRUE to merge on grid
  // after re-running the jobs in SetRunMode("terminate") 
  // (see below) mode, set SetMergeViaJDL(kFALSE) 
  // to collect final results
  alienHandler->SetOneStageMerging(kFALSE); // One can use this to force
  alienHandler->SetMaxMergeStages(2);
  alienHandler->SetMergeViaJDL(kTRUE);

  // define the output folders
  alienHandler->SetGridWorkingDir("myTestDir");
  alienHandler->SetGridOutputDir("myTestOutDir");
  alienHandler->SetRunMode("test");
  alienHandler->SetNtestFiles(1);

  alienHandler ->SetDefaultOutputs();
  alienHandler->SetOutputToRunNo();
  alienHandler->SetAnalysisMacro("vnRun2.C");
  alienHandler->SetSplitMaxInputFileNumber(10);
  alienHandler->SetMasterResubmitThreshold(90);
  alienHandler->SetTTL(30000);
  alienHandler->SetInputFormat("xml-single");
  alienHandler->SetJDLName("vnRun2.jdl");
  alienHandler->SetExecutable("vnRun2.sh");
  alienHandler->SetPrice(1);    
  alienHandler->SetSplitMode("se");
  //alienHandler->SetExecutableCommand("aliroot -b -q");

  // connect the alien alienHandler to the manager
  mgr->SetGridHandler(alienHandler);
 
  //Add this here: run before your task, but after definition of manager and input handler
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C");
  AliMultSelectionTask * task = AddTaskMultSelection(kFALSE); // user mode:
  
  AliAnalysisTaskChargedFlow* taskFlow1 = AddTaskChargedFlow(128, 0.8, 10.0, 0.2, 5.0, 2, 1, 70, 0, 1, kTRUE, kTRUE, 14, 0, kFALSE, "TPConly", kFALSE);
  taskFlow1->SelectCollisionCandidates(AliVEvent::kINT7); //..V0 offline selection
   

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();


  // and launch the analysis
  mgr->StartAnalysis("grid");





}
