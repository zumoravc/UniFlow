#!/bin/bash

outPath=~/NBI/Codes/results/V0s/5/plusplus/
cd ${outPath}

mkdir -pv plots
mkdir -pv plots/InvMassK0s
mkdir -pv plots/InvMassLambda
mkdir -pv plots/CummMassK0s
mkdir -pv plots/CummMassLambda
mkdir -pv plots/FlowMassK0s
mkdir -pv plots/FlowMassLambda

root -l -b -q ~/NBI/Codes/macros/ProcessV0s.C\(\"${outPath}/merge/AnalysisResults_merged.root\",\"${outPath}/plots\"\)