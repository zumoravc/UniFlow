// shortcut for running single macro

void runSystPlot()
{
  gROOT->Macro("Load.C");

  TString sTaskTag;

  // sTaskTag = "FB768";
  // sTaskTag = "PID3sigma";
  // sTaskTag = "PVz8";
  // sTaskTag = "TPCcls90";
  // sTaskTag = "V0sCPA099";
  // sTaskTag = "V0sCrossFind1";
  // sTaskTag = "V0sDaugDCA3";
  // sTaskTag = "V0sDaugPt02";
  // sTaskTag = "V0sDecRad10";
  // sTaskTag = "V0sFinderOn";
  sTaskTag = "V0sPVDCA3";

  gROOT->Macro(Form("RunMixedSyst.C(\"%s\")",sTaskTag.Data()));
  gROOT->Macro(Form("PlotFits.C(\"../results/nlf/systematics/Lambda/%s/\")",sTaskTag.Data()));
}
