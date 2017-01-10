#!/bin/bash

# Cleaning directory from task submission generated files 
cd ~/NBI/Flow/flow/

rm ./*.xml
rm ./FlowPID*
rm ./myAnalysis.C
rm ./event_stat.root

# compilation files
rm ./*.d
rm ./*.so


