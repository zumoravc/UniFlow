#!/bin/bash

# ===========================================================
# Script copying code to path tag
# Arguments:
#		1 - tag
#   2 - runlist
# ===========================================================

# checking the parameters
if [ $# -ne 2 ]; then
	echo "Wrong number of parameters: $# (2 needed). Exit!"
	exit
fi

if [ "$1" = "" ] || [ "$2" = "" ]; then
	echo "Empty parameter given. Exit!"
	exit
fi

if [ ! -f "/Users/vpacik/NBI/Flow/runLists/$2.runlist" ]; then
	echo "Runlist file (arg. 2) does not exists. Exit!"
	exit
fi

# passing checks

path=~/NBI/Flow/results/${1}
code=~/NBI/Flow/flow/

#echo "Argument: \"$1\"" 
#echo "Full path:${path}"

# copying files into $path
echo "Copying files from ${code} to ${path}:"
mkdir -pv ${path}

# copy code files
cp -v ${code}/AliAnalysisTaskFlowPID.cxx ${path}
cp -v ${code}/AliAnalysisTaskFlowPID.h ${path}
cp -v ${code}/AddTaskFlowPID.C ${path}
cp -v ${code}/runPP.C ${path}
cp -v ${code}/runPPb.C ${path}

# generate maintenance scripts
. ${code}/scripts/generateMergeOnGrid.sh ${path} ${2} ${1}

#cp -v ${code}/scripts/CleanAfterMerge.sh ${path}
#cp -v ${code}/scripts/SubmitMergeOnGrid.sh ${path}
#cp -v ${code}/scripts/SubmitMergeLocal.sh ${path}

## change directory to $path
#echo "Change directory to ${path}"


