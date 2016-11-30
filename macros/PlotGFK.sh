#!/bin/bash

inputFile=~/NBI/Flow/results/5-GFK-PbPb-MultSelectNoRejection/Flow.root
outPath=~/NBI/Flow/results/5-GFK-PbPb-MultSelectNoRejection/comp-try2/

mkdir -pv ${outPath}
mkdir -pv ${outPath}/Ref
mkdir -pv ${outPath}/Tracks
mkdir -pv ${outPath}/Tracks/Ratios
mkdir -pv ${outPath}/Pions
mkdir -pv ${outPath}/Kaons
mkdir -pv ${outPath}/Protons

root -b -l  ~/NBI/Flow/macros/PlotGFK.C\(\"${outPath}\",\"${inputFile}\",1,1,1\)
