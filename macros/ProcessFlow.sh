#!/bin/bash

outDir="~/NBI/Codes/results/TPConly"

macroPath="~/NBI/Codes/macros"

mkdir -pv ${outDir}
mkdir -pv ${outDir}/CompYouRef
mkdir -pv ${outDir}/CompKatarinaDiff

alienv load AliPhysics/latest

root -l -b -q ${macroPath}/ProcessFlow.C
