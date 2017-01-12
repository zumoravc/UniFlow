#!/bin/bash

# ===========================================================
# Script copying code to path tag
# Arguments:
#   1 - "-system" (-pp / -ppb -pbpb)
#		2 - tag (witout the runlist)
#   3 - runlist
# ===========================================================

# checking the parameters
if [ $# -ne 3 ]; then
	echo "Wrong number of parameters: $# (3 needed). Exit!"
	exit
fi

if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ]; then
	echo "Empty parameter given. Exit!"
	exit
fi

case ${1} in
	-pp|-PP )
	echo "Analysis on PP"
	system="pp"
	;;
	-pPb|-PPB|-ppb|-PPb )
	echo "Analysis on pPb"
	system="pPb"
	;;
	-pbpb|-PbPb|-PBPB )
	echo "Analysis on PbPb"
	system="PbPb"
	;;
	* )
	echo "System not selected; Exit!"
	exit
	;;
esac

if [ ! -f "/Users/vpacik/NBI/Flow/runLists/$3.runlist" ]; then
	echo "Runlist file '${3}.runlist' (arg. 3) does not exists. Exit!"
	exit
fi
tag=${system}-${3}-${2}
path=~/NBI/Flow/results/${tag}
code=~/NBI/Flow/flow/

if [ -d "${path}" ]; then
	# check if directory given by tag does not exists
	echo "Output directory '${path}' already exists!"
	exit
fi

# passing checks

#echo "Argument: \"$1\""
#echo "Full path:${path}"

# copying files into $path
echo "Copying files from ${code} to ${path}:"
mkdir -pv ${path}

# copy code files
cp -v ${code}/AliAnalysisTaskFlowPID.cxx ${path}
cp -v ${code}/AliAnalysisTaskFlowPID.h ${path}
cp -v ${code}/AddTaskFlowPID.C ${path}
#cp -v ${code}/runMC.C ${path}

if [ "${system}" = "pp" ]; then
	cp -v ${code}/runPP.C ${path}
fi

if [ "${system}" = "pPb" ]; then
	cp -v ${code}/runPPb.C ${path}
fi

# generate maintenance scripts
. ${code}/scripts/generateMergeOnGrid.sh ${path} ${3} ${tag}
. ${code}/scripts/generateMergeLocal.sh ${path}

#cp -v ${code}/scripts/CleanAfterMerge.sh ${path}

## change directory to $path
#echo "Change directory to ${path}"
exit
