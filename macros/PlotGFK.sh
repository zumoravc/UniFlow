#!/bin/bash

inputFile=~/NBI/Flow/results/1-GFK-PbPb-YouBinning/Flow-YouBins.root
outPath=~/NBI/Flow/results/1-GFK-PbPb-YouBinning/comp

mkdir -pv ${outPath}
mkdir -pv ${outPath}/Ref
mkdir -pv ${outPath}/Tracks
mkdir -pv ${outPath}/Tracks/Ratios
mkdir -pv ${outPath}/Pions
mkdir -pv ${outPath}/Kaons
mkdir -pv ${outPath}/Protons

root -l -q ~/NBI/Flow/macros/PlotGFK.C\(\"${outPath}\",\"${inputFile}\"\)
