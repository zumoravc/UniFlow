#!/bin/bash

# ===========================================================
# Script for generating SubmitMergeLocal.sh to given folder
# Arguments:
#		1 - output folder where generate
# ===========================================================

echo "##### Generating mergeLocal.sh #####"

# checking parameters
if [ $# -ne 1 ]; then
	echo "Wrong number of parameters: $# (1 needed). Exit!"
	exit
fi

if [ "$1" = "" ] ; then
	echo "Empty parameter given. Exit!"
	exit
fi

output=$1

if [ ! -d "${output}" ]; then
	#path folder does not exist
  echo "Output folder (arg. 1) does not exists. Exit!"
	exit
fi

echo "--- Listing arguments ---"
echo "output:\"${output}\""
echo "-------------------------"

# writing following script into $file
cat > ${output}/mergeLocal.sh <<EOT
#!/bin/bash
outPath=${output}

cd \${outPath}/merge

echo "Checking if local merging result file exits (file: ${outPath}/merge/AnalysisResults.root)"

if [ ! -f ./AnalysisResults.root ]; then
	echo "File NOT found! Merging locally."
	root -l -b -q ~/NBI/Flow/macros/mergeOutput.C
fi

echo "Checking if file exists now"

if [ ! -f ./AnalysisResults.root ]; then
	echo "File NOT found! Something went wrong. EXITING"
	exit
fi

echo "File exists. Cleaning!"

cd ../
mv -v ./merge/AnalysisResults.root ./

if [ -f ./AnalysisResults.root ]; then
	rm -rf ./merge
	rm ./*.xml
	echo "Cleaning succesful!"
fi
EOT

if [ -f ${output}/mergeLocal.sh ]; then
	chmod +x ${output}/mergeLocal.sh
	echo "File generated!"
else
	echo "File NOT generated!"
	exit
fi
