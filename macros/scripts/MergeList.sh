#!/bin/bash

runList="10eAOD135full.txt"
prefix="000"
output="LHC10e/AnalysisResults10e.root"

input=""
for i in `cat $runList`
do
  path="LHC10e/$prefix$i/AnalysisResults.root"
  input="$input $path"
done

echo $input

hadd $output $input
