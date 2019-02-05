#!/bin/bash

# script for clearing redundant file after running on GRID & after succesfull local merging

rm -rf ./tempChecksum_*
rm ./*.pcm
rm ./*.xml
rm ./*_cxx.d
rm ./*_cxx.so
rm ./*ACLiC_cxx.*
rm ./UniFlow*
rm ./myAnalysis.C
rm ./FlowPID*
rm ./gdb-backtrace*
rm ./std*
rm ./outputs_valid

echo "Cleaning done!"
