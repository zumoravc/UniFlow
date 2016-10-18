#!/bin/bash
inPath=~/NBI/Flow/results/V0s/11-filtering/merge
outPath=~/NBI/Flow/results/V0s/11-filtering
tag=JHEP

# =================================================

mkdir -pv ${outPath}
cd ${outPath}

for gap in 08 #10 #02 04 06 08 10
do
	echo " === Processing Gap${gap} === "
	plotsDir=plots_${tag}_Gap${gap}
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
	mkdir -pv ${plotsDir}/finalFlow/


	#root -l -b -q ~/NBI/Flow/macros/ProcessV0s.C\(\"${inPath}/AnalysisResults.root\",\"${outPath}/${plotsDir}\",\"flowPID_${tag}\",\"Gap${gap}\",\"png\"\)
	root -l -b -q ~/NBI/Flow/macros/V0sExtractFlow.C\(\"${outPath}/${plotsDir}/V0sFlow_Gap${gap}.root\",\"${outPath}/${plotsDir}\",\"flowPID_${tag}\",\"Gap${gap}\",\"png\",kTRUE\)
done
