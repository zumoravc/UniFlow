#!/bin/bash
outPath=~/NBI/Flow/results/V0s/QApid/

cd ${outPath}/merge

root -l -b -q ~/NBI/Flow/macros/mergeOutput.C
