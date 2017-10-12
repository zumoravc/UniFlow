#!/bin/bash

# script for clearing redundant file after running on GRID & after succesfull local merging

rm ./*.xml
rm ./*_cxx.d
rm ./*_cxx.so
rm ./UniFlow*
rm ./myAnalysis.C

# checking if AnalysisResults.root in "root dir" exits, if so, deleting the per-run files
if [ -f ./AnalysisResults.root ]; then
	echo "File './AnalysisResults.root' found! Removing the 'merge' dir!"
  rm -rf ./merge/
fi

rm ./mergeOnGrid.sh
rm ./mergeLocal.sh
rm ./cleanAfterMerge.sh

echo "Cleaning done!"
