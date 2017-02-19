#!/bin/bash

# ===========================================================
# Script copying code to path tag
# Arguments:
#   1 - system (pp / ppb pbpb)
#   2 - period (only full)
#		3 - tag (witout the runlist)
# ===========================================================

# checking the parameters
if [ $# -ne 3 ]; then
	echo "Wrong number of parameters: $# (3 needed). Exit!"
	echo "Arguments: system period tag !"
	exit
fi

if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ]; then
	echo "Empty parameter given. Exit!"
	exit
fi

case ${1} in
	pp|PP )
	echo "Analysis on PP"
	system="pp"
	;;
	pPb|PPB|ppb|PPb )
	echo "Analysis on pPb"
	system="pPb"
	;;
	pPbRun1|PPBRun1|ppbrun1|PPbRun1 )
	echo "Analysis on pPb (Run1)"
	system="pPbRun1"
	;;
	pPbRun2|PPBRun2|ppbrun2|PPbRun2 )
	echo "Analysis on pPb (Run2)"
	system="pPbRun2"
	;;
	pbpb|PbPb|PBPB )
	echo "Analysis on PbPb"
	system="PbPb"
	;;
	* )
	echo "System not selected; Exit!"
	exit
	;;
esac

if [ ! -f "/Users/vpacik/NBI/Flow/runLists/$2.runlist" ]; then
	echo "Runlist file '${2}.runlist' (arg. 3) does not exists. Exit!"
	exit
fi

period=${2}
desc=${3}

tag=${desc}-${system}-${period}
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
	#cp -v ${code}/runPP.C ${path}
	. ${code}/scripts/generateRunPP.sh ${path} ${tag} ${period}
fi

if [ "${system}" = "pPb" ]; then
	cp -v ${code}/runPPb.C ${path}
fi

if [ "${system}" = "pPbRun1" ]; then
	cp -v ${code}/runPPbRun1.C ${path}
fi

if [ "${system}" = "pPbRun2" ]; then
	#cp -v ${code}/runPPbRun2.C ${path}
	. ${code}/scripts/generateRunPPbRun2.sh ${path} ${tag} ${period}
fi

# generate maintenance / runnning scripts
. ${code}/scripts/generateMergeOnGrid.sh ${path} ${period} ${tag}
. ${code}/scripts/generateMergeLocal.sh ${path}

#if [ "${system}" = "pp" ]; then
#	. ${code}/scripts/generateRunPP.sh ${path} ${tag} ${period}
#fi

#cp -v ${code}/scripts/CleanAfterMerge.sh ${path}

## change directory to $path
#echo "Change directory to ${path}"
exit
