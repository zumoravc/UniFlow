#!/bin/bash

inputFile=~/NBI/Flow/results/0-GFK-PbPb-full-10samples/Flow.root
outPath=~/NBI/Flow/results/0-GFK-PbPb-full-10samples/comp2

mkdir -pv ${outPath}
mkdir -pv ${outPath}/Ref
mkdir -pv ${outPath}/Tracks
mkdir -pv ${outPath}/Tracks/Ratios
mkdir -pv ${outPath}/Pions
mkdir -pv ${outPath}/Kaons
mkdir -pv ${outPath}/Protons

root -b -l -q ~/NBI/Flow/macros/PlotGFK.C\(\"${outPath}\",\"${inputFile}\",1,0,1\)
