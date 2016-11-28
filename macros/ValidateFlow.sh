#!/bin/bash

#inputFile=~/NBI/Flow/results/7-GFK-PbPb-OldCent2/Flow_PileOFFPeriodOFF.root
inputFile=~/NBI/Flow/results/8-GFK-PbPb-UltCentTest/Flow_FB768_New_PileON_PeriodOFF.root
outPath=~/NBI/Flow/results/8-GFK-PbPb-UltCentTest/validate_ver3_New_ON_OFF

mkdir -pv ${outPath}
mkdir -pv ${outPath}/Ref
mkdir -pv ${outPath}/Tracks
mkdir -pv ${outPath}/Tracks/Ratios
mkdir -pv ${outPath}/Pions
mkdir -pv ${outPath}/Kaons
mkdir -pv ${outPath}/Protons

root -b -l -q ~/NBI/Flow/macros/ValidateFlow.C\(\"${outPath}\",\"${inputFile}\",1,1,1\)
