#!/bin/bash
inPath=~/NBI/Codes/results/V0s/8/plusplus_part1/merge
outPath=~/NBI/Codes/results/V0s/8/plusplus_part1

mkdir -pv ${outPath}
cd ${outPath}

mkdir -pv plots
# ProcessV0s.C requirements
mkdir -pv plots/InvMassK0s
mkdir -pv plots/InvMassLambda
mkdir -pv plots/CummMassK0s
mkdir -pv plots/CummMassLambda
mkdir -pv plots/FlowMassK0s
mkdir -pv plots/FlowMassLambda
mkdir -pv plots/compInvMass/
mkdir -pv plots/compFlowMass/

# V0sExtractFlow.C requirements
mkdir -pv plots/fitK0s/
mkdir -pv plots/fitLambda/

#root -l -b -q ~/NBI/Codes/macros/ProcessV0s.C\(\"${inPath}/AnalysisResults_merged.root\",\"${outPath}/plots\",\"flowPID_JHEP\",\"Gap00\",\"png\"\)
root -l -b -q ~/NBI/Codes/macros/V0sExtractFlow.C\(\"${outPath}/plots/V0sFlow.root\",\"${outPath}/plots\",\"flowPID_JHEP\",\"Gap00\",\"png\"\)
