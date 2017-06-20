#!/bin/bash

outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst/tpcCls_2/tpcCls80"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst/NUA_cor/CENT_woSDD_16q/bgFunc"

echo "Prepring folder structure for UniFlow in '${outputPath}'"

mkdir -v $outputPath

# preparing folder for binning suggestion
mkdir -v ${outputPath}/suggestBins
mkdir -v ${outputPath}/plots

# preparing folder for slices
# mkdir -pv ${outputPath}/fits/Refs
mkdir -pv ${outputPath}/slices/Charged
mkdir -pv ${outputPath}/slices/Pion
mkdir -pv ${outputPath}/slices/Kaon
mkdir -pv ${outputPath}/slices/Proton
mkdir -pv ${outputPath}/slices/Phi
mkdir -pv ${outputPath}/slices/K0s
mkdir -pv ${outputPath}/slices/Lambda

# preparing folder for fits
mkdir -pv ${outputPath}/fits/Phi
mkdir -pv ${outputPath}/fits/K0s
mkdir -pv ${outputPath}/fits/Lambda

echo "Preparation done"
