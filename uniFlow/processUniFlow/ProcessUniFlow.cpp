/* ProcessUniFlow class
 *
 * Class implemented for processing results of AliAnalysisTaskUniFlow task.
 *
 * Author: Vojtech Pacik (vojtech.pacik@cern.ch), NBI, 2017
 */

#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TString.h"

class ProcessUniFlow
{
  public:
                ProcessUniFlow();
                ~ProcessUniFlow();

    void        SetInputFilePath(const char* path) { fsInputFilePath = path; }
    void        SetInputFileName(const char* name) { fsInputFileName = name; }
    void        SetOutputFilePath(const char* path) { fsOutputFilePath = path; }
    void        SetOutputFileName(const char* name) { fsOutputFileName = name; }
    void        SetDebug(Bool_t debug = kTRUE) { fbDebug = debug; }
    void        Run(); // running the task (main body of the class)

  protected:

  private:
    Bool_t      Initialize(); // initialization task

    // printing output methods
    void        Fatal(TString sMsg, TString sMethod = ""); // printf the msg as error
    void        Error(TString sMsg, TString sMethod = ""); // printf the msg as error
    void        Warning(TString sMsg, TString sMethod = ""); // printf the msg as warning
    void        Info(TString sMsg, TString sMethod = ""); // printf the msg as info
    void        Debug(TString sMsg, TString sMethod = ""); // printf the msg as info

    Bool_t      fbInit; // flag for initialization status
    Bool_t      fbDebug; // flag for debugging : if kTRUE Debug() messages are displayed
    TFile*      ffInputFile; //! input file container
    TFile*      ffOutputFile; //! input file container
    TString     fsInputFilePath; // path to the input folder with input file
    TString     fsInputFileName; // name of input file
    TString     fsOutputFilePath; // path to the ouput folder
    TString     fsOutputFileName; // name of output file
};
//_____________________________________________________________________________
ProcessUniFlow::ProcessUniFlow() :
  fbDebug(kFALSE),
  fbInit(kFALSE),
  ffInputFile(0x0),
  ffOutputFile(0x0)
{
  // default constructor
  fsInputFilePath = TString("");
  fsInputFileName = TString("AnalysisResults.root");
  fsOutputFilePath = TString("");
  fsOutputFileName = TString("UniFlow.root");
}
//_____________________________________________________________________________
ProcessUniFlow::~ProcessUniFlow()
{
  // default destructor
  if(ffInputFile) delete ffInputFile;
  if(ffOutputFile) delete ffOutputFile;
}
//_____________________________________________________________________________
void ProcessUniFlow::Run()
{
  // main body of the class
  if(!Initialize()) { Fatal("Task not initialized","Run"); return; }

  Debug("Initialized");

  return;
}
//_____________________________________________________________________________
Bool_t ProcessUniFlow::Initialize()
{
  // initialization of all necessery prerequisits
  Info("Initializating task","Initialize");
  fbInit = kFALSE;

  // opening input file
  ffInputFile = new TFile(Form("%s/%s",fsInputFilePath.Data(),fsInputFileName.Data()),"READ");
  if(!ffInputFile || !ffInputFile->IsOpen())
  {
    Fatal(Form("Input file %s/%s not open",fsInputFilePath.Data(),fsInputFileName.Data()),"Initialize");
    return fbInit;
  }

  // opening output file
  ffOutputFile = new TFile(Form("%s/%s",fsOutputFilePath.Data(),fsOutputFileName.Data()),"RECREATE");
  if(!ffOutputFile || !ffOutputFile->IsOpen())
  {
    Fatal(Form("Output file %s/%s not open",fsOutputFilePath.Data(),fsOutputFileName.Data()),"Initialize");
    return fbInit;
  }


  // initialization succesfull
  fbInit = kTRUE;
  Info("Initialization succesfull","Initialize");
  return fbInit;
}
//_____________________________________________________________________________
void ProcessUniFlow::Fatal(TString sMsg, TString sMethod)
{
	printf("Fatal::%s  %s. Terminating!\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Error(TString sMsg, TString sMethod)
{
	printf("Error::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Info(TString sMsg, TString sMethod)
{
	printf("Info::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Warning(TString sMsg, TString sMethod)
{
	printf("Warning::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
void ProcessUniFlow::Debug(TString sMsg, TString sMethod)
{
	if(fbDebug) printf("Debug::%s  %s\n", sMethod.Data(), sMsg.Data());
}
//_____________________________________________________________________________
