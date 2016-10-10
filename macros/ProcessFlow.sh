#!/bin/bash

inDir=~/NBI/Codes/results/V0s/8/merge/plusplus
outDir=~/NBI/Codes/results/V0s/8/JHEP

macroPath=~/NBI/Codes/macros

mkdir -pv ${outDir}
mkdir -pv ${outDir}/CompYouRef
mkdir -pv ${outDir}/CompYouDiff
mkdir -pv ${outDir}/CompKatarinaDiff

#alienv enter AliPhysics/latest

root -l -b -q ${macroPath}/ProcessFlow.C\(\"${inDir}/AnalysisResults.root\",\"${outDir}\"\)
