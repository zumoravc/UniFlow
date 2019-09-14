#!/bin/bash
tag=pp-2016l-FB32

path=~/NBI/Flow/results

cd ${path}/${tag}/

if [ ! -f ./merge/AnalysisResults.root ]; then
	echo "File not found"
	exit
fi

mv -v ./merge/AnalysisResults.root ./
rm -rf ./merge
