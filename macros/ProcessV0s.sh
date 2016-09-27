#!/bin/bash
inPath=~/NBI/Codes/flow
outPath=~/NBI/Codes/results/V0s/test_tag/plusplus_JDL
mkdir -pv ${outPath}
cd ${outPath}

mkdir -pv plots
mkdir -pv plots/InvMassK0s
mkdir -pv plots/InvMassLambda
mkdir -pv plots/CummMassK0s
mkdir -pv plots/CummMassLambda
mkdir -pv plots/FlowMassK0s
mkdir -pv plots/FlowMassLambda
mkdir -pv plots/compInvMass/
mkdir -pv plots/compFlowMass/

root -l -b -q ~/NBI/Codes/macros/ProcessV0s.C\(\"${inPath}/AnalysisResults.root\",\"${outPath}/plots\",\"flowPID_JHEP\"\)