#!/bin/bash

# ===========================================================
# Script for download merged output (run-by-run) from Grid
# ===========================================================

tag=pp-16l-run2
runList=/Users/vpacik/NBI/Flow/runLists/lhc16l_v2.runlist

mkdir -pv merge
cd merge

for i in $(cat ${runList})
do
	mkdir -pv merge_${i}
	cd merge_${i}
	alien_cp alien:/alice/cern.ch/user/v/vpacik/${tag}/output/000${i}/AnalysisResults.root ./
	cd ../
done
