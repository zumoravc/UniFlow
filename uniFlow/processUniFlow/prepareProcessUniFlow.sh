#!/bin/bash

# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTestKaon"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTestKaon"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/filterBit_withNUA/UniFlow_fb768_kaons"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/tpcCls_merged/uniflow_tpcCls80"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/v0s/UniFlow_DecayRadius"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/Vtx_z/UniFlow_vtx8"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/PbPblike"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/uniflow_NoArmenteros"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/kaonDist/Alex"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/loosePhi_tightV0s/kaon"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/fit_pol/fit_pol"
outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/NUA_cor/merged_16qt_noFASTt/fitTest_afterFinal"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/pid_merged/"
# outputPath="/Users/vpacik/NBI/Flow/results/uniFlow_syst_run2/Vtx_z_merged/UniFlow_vtx9"

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
