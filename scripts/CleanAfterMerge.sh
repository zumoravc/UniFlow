#!/bin/bash
tag=11-GFK-PbPb-R1-BayesTest

path=~/NBI/Flow/results

cd ${path}/${tag}/
mv ./merge/AnalysisResults.root ./
rm -rfv ./merge