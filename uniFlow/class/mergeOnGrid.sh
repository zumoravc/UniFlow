#!/bin/bash

# ===========================================================
# Script for local merging & download output on Grid
# ===========================================================

tag=pp-16l-run2
runList=/Users/vpacik/NBI/Flow/runLists/lhc16l_v2.runlist

mkdir -pv merge
cd merge

for i in $(cat ${runList})
do
	mkdir -pv merge_${i}
	cd merge_${i}
	path="/alice/cern.ch/user/v/vpacik/${tag}/output/000${i}/"

	root -l -b -q ~/NBI/Flow/macros/mergeOutputOnGrid.C\(\"${path}\"\) &
	# alien_cp alien:${path}AnalysisResults.root ./

	cd ../
done
