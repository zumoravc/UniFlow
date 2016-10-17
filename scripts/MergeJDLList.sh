#!/bin/bash

runList="10eAOD135test.txt"
prefix="000"
iter=1

for i in `cat $runList`
do
  if [ $iter -eq 6 ]
  then
    echo "=== Submitted $iter iterations. Sleep for 2 minutes ==="
    sleep 120
    iter=1
  fi
  iter=`expr $iter + 1`

  echo "=== Merging job $i ==="
  path="LHC10e/$prefix$i"
  #mkdir -pv $path
  #cp -r ./template/* ./$path/
  cd $path/
  root -l -b -q "runEMCalJetGrid.C(kTRUE,kTRUE,\"terminate\",\"$i\")" > stdMerge.log 2>&1 &
  cd ../../
done
