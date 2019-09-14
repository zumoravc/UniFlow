#!/bin/bash
outPath=~/NBI/Flow/results/V0s/QA/

mkdir -pv ${outPath}
cd ${outPath}
mkdir -pv merge
cd merge

for i in $(cat ~/NBI/Flow/runLists/lhc10h_plusplus.runlist)
do
	mkdir -pv merge_${i}
	cd merge_${i}
	path="/alice/cern.ch/user/v/vpacik/V0s/QApid/outFlow/000${i}/"

	#echo $path

	root -l -b -q ~/NBI/Flow/macros/mergeOutputOnGrid.C\(\"${path}\"\) &
	cd ../

done
