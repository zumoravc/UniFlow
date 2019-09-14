#!/bin/bash

runList="10eAOD135full.txt"
prefix="000"
iter=1

for i in `cat $runList`
do
  if [ $iter -eq 6 ]
  then
    echo "=== Submitted $iter iterations. Sleep for 4 minutes ==="
    sleep 240
    iter=1
  fi
  iter=`expr $iter + 1`

  echo "=== Submitting $prefix$i ==="
  path="LHC10e/$prefix$i"
  mkdir -pv $path
  cp -r ./template/* ./$path/
  cd $path/
  root -l -b -q "runEMCalJetGrid.C(kTRUE,kTRUE,\"full\",\"$i\")" > stdSub.log 2>&1 &
  cd ../../
done
