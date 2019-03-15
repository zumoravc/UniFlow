// shortcut for running single macro

void runSyst()
{
  gROOT->Macro("Load.C");
  gROOT->Macro("RunMixedSyst.C");
}
