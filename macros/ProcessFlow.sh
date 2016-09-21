#!/bin/bash

outDir="/home/vojtech/NBI/Codes/results/TPConly"

macroPath="/home/vojtech/NBI/Codes/macros"

mkdir -pv ${outDir}
mkdir -pv ${outDir}/CompYouRef
mkdir -pv ${outDir}/CompYouDiff
mkdir -pv ${outDir}/CompKatarinaDiff

#alienv enter AliPhysics/latest

root -l -b -q ${macroPath}/ProcessFlow.C
