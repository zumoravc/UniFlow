#!/bin/bash

#inputFile=~/NBI/Flow/results/7-GFK-PbPb-OldCent2/Flow_PileOFFPeriodOFF.root
inputFile=~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/Flow_FB768_Old_PileON_PeriodON.root
outPath=~/NBI/Flow/results/10-GFK-PbPb-R1-noSampling/validate_Old_ON_ON

mkdir -pv ${outPath}
mkdir -pv ${outPath}/Ref
mkdir -pv ${outPath}/Tracks
mkdir -pv ${outPath}/Tracks/Ratios
mkdir -pv ${outPath}/Pions
mkdir -pv ${outPath}/Kaons
mkdir -pv ${outPath}/Protons

root -b -l -q ~/NBI/Flow/macros/ValidateFlow.C\(\"${outPath}\",\"${inputFile}\",1,1,1\)
