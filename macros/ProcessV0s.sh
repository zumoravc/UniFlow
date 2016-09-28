#!/bin/bash
inPath=~/NBI/Codes/results/V0s/7/plusplus
outPath=~/NBI/Codes/results/V0s/7/plusplus

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

#root -l -b -q ~/NBI/Codes/macros/ProcessV0s.C\(\"${inPath}/AnalysisResults_merged.root\",\"${outPath}/plots\",\"flowPID_JHEP\"\)
root -l -b -q ~/NBI/Codes/macros/ProcessV0s.C\(\"${inPath}/AnalysisResults_merged.root\",\"${outPath}/plots\",\"flowPID_lose\"\)