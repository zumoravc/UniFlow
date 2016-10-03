#!/bin/bash
outPath=~/NBI/Codes/results/V0s/8/

mkdir ${outPath}
cd ${outPath}
mkdir -pv merge
cd merge

for i in $(cat ~/NBI/Codes/runLists/lhc10h_minmin_part1.runlist)
do
	mkdir -pv merge_${i}
	cd merge_${i}
	path="/alice/cern.ch/user/v/vpacik/V0s/8_minusminus_part1/outFlow/000${i}/"
	
	#echo $path

	root -l -b -q ~/NBI/Codes/macros/mergeOutputOnGrid.C\(\"${path}\"\) &
	cd ../

done
