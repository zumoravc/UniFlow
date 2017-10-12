#!/bin/bash

# ===========================================================
# Script for local merging & download output on Grid
# ===========================================================

tag=K0sComp_0
runList=/Users/vpacik/NBI/Flow/runLists/temp.runlist

mkdir -pv merge
cd merge

for i in $(cat ${runList})
do
	mkdir -pv merge_${i}
	cd merge_${i}
	path="/alice/cern.ch/user/v/vpacik/${tag}/output/00${i}/"

	root -l -b -q ~/NBI/Flow/macros/mergeOutputOnGrid.C\(\"${path}\"\) &
	cd ../
done
