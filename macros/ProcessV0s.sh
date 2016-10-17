#!/bin/bash
inPath=~/NBI/Flow/results/V0s/10/merge/
outPath=~/NBI/Flow/results/V0s/10
tag=JHEP

# =================================================

mkdir -pv ${outPath}
cd ${outPath}

plotsDir=plots_${tag}
mkdir -pv ${plotsDir}
# ProcessV0s.C requirements
mkdir -pv ${plotsDir}/InvMassK0s
mkdir -pv ${plotsDir}/InvMassLambda
mkdir -pv ${plotsDir}/CummMassK0s
mkdir -pv ${plotsDir}/CummMassLambda
mkdir -pv ${plotsDir}/FlowMassK0s
mkdir -pv ${plotsDir}/FlowMassLambda
mkdir -pv ${plotsDir}/compInvMass/
mkdir -pv ${plotsDir}/compFlowMass/

# V0sExtractFlow.C requirements
mkdir -pv ${plotsDir}/fitK0s/
mkdir -pv ${plotsDir}/fitLambda/

#root -l -b -q ~/NBI/Flow/macros/ProcessV0s.C\(\"${inPath}/AnalysisResults.root\",\"${outPath}/${plotsDir}\",\"flowPID_${tag}\",\"Gap00\",\"png\"\)
root -l -b -q ~/NBI/Flow/macros/V0sExtractFlow.C\(\"${outPath}/${plotsDir}/V0sFlow.root\",\"${outPath}/${plotsDir}\",\"flowPID_${tag}\",\"Gap00\",\"png\"\)