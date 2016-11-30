#!/bin/bash
tag=FB768
inPath=~/NBI/Flow/results/0-GFK-PbPb-full/merge
outPath=~/NBI/Flow/results/0-GFK-PbPb-full/${tag}

# =================================================

mkdir -pv ${outPath}
cd ${outPath}

for gap in 00 # 00 10 #10 #02 04 06 08 10
do
	echo " === Processing Gap${gap} === "
	plotsDir=plots_Gap${gap}
	mkdir -pv ${plotsDir}
	# ProcessV0s.C requirements
	mkdir -pv ${plotsDir}/InvMassK0s
	mkdir -pv ${plotsDir}/InvMassLambda
	mkdir -pv ${plotsDir}/CumMassK0s
	mkdir -pv ${plotsDir}/CumMassLambda
	mkdir -pv ${plotsDir}/FlowMassK0s
	mkdir -pv ${plotsDir}/FlowMassLambda
	mkdir -pv ${plotsDir}/RefFlow
	
	# V0sExtractFlow.C requirements
	mkdir -pv ${plotsDir}/fitK0s/
	mkdir -pv ${plotsDir}/fitLambda/
	mkdir -pv ${plotsDir}/PtFlow/


	root -l -b -q ~/NBI/Flow/macros/ProcessV0s.C\(\"${inPath}/AnalysisResults.root\",\"${outPath}/${plotsDir}\",\"flowPID_${tag}\",\"Gap${gap}\",\"png\"\)
	#root -l -b -q ~/NBI/Flow/macros/ProcessV0sSampled.C\(\"${inPath}/AnalysisResults.root\",\"${outPath}/${plotsDir}\",\"flowPID_${tag}\",\"Gap${gap}\",\"png\"\)
	
	#root -l -b -q ~/NBI/Flow/macros/V0sExtractFlow.C\(\"${outPath}/MassDist_V0s_Gap${gap}.root\",\"${outPath}/${plotsDir}\",\"flowPID_${tag}\",\"Gap${gap}\",\"png\"\)
done