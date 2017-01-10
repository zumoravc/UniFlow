#!/bin/bash

# ===========================================================
# Script for generating SubmitMergeOnGrid.sh to given folder
# Arguments:
#		1 - output folder where generate
#   2 - runlist
#   3 - tag
# ===========================================================

echo "##### Generating mergeOnGrid.sh #####"

# checking parameters
if [ $# -ne 3 ]; then
	echo "Wrong number of parameters: $# (3 needed). Exit!"
	exit
fi

if [ "$1" = "" ] || [ "$2" = "" ] || [ "$3" = "" ] ; then
	echo "Empty parameter given. Exit!"
	exit
fi

output=$1
run=~/NBI/Flow/runLists/$2.runlist
task=$3

if [ ! -d "${output}" ]; then 
	#path folder does not exist
  echo "Output folder (arg. 1) does not exists. Exit!"
	exit
fi

if [ ! -f "${run}" ]; then 
	# runlist file does not exist
  echo "Runlist file (arg. 2) does not exists. Exit!"
	exit
fi
# arguments passed

echo "--- Listing arguments ---"
echo "output:\"${output}\""
echo "runlist:\"${run}\""
echo "tag:\"${task}\""
echo "-------------------------"


# writing following script into $file
cat > ${output}/mergeOnGrid.sh <<EOT
#!/bin/bash

# ===========================================================
# Script for local merging & download output on Grid 
# ===========================================================

tag=${task}
runList=${run}

mkdir -pv merge
cd merge

for i in \$(cat \${runList})
do
	mkdir -pv merge_\${i}
	cd merge_\${i}
	path="/alice/cern.ch/user/v/vpacik/\${tag}/outFlow/000\${i}/"
	
	root -l -b -q ~/NBI/Flow/macros/mergeOutputOnGrid.C\(\"\${path}\"\) &
	cd ../
done
EOT

if [ -f ${output}/mergeOnGrid.sh ]; then
	chmod +x ${output}/mergeOnGrid.sh
	echo "File generated!"
else
	echo "File NOT generated!"
	exit
fi
