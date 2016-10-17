#!/bin/bash
runList="10eAOD135test.txt"
prefix="000"
sleepSec=20
echo "### Submitting runs from $runList ###"

for i in `cat $runList`
do
  echo "=== Submitting $prefix$i ==="
  path="test/$prefix$i"
  mkdir -pv $path
  cp -r ./template/* ./$path/
  cd $path/
  root -l -b -q "runEMCalJetGrid.C(kTRUE,kTRUE,\"full\",\"$i\")" > stdSub.log 2>&1 &
  cd ../../
  sleep 20
done
echo "### Submitting DONE ###"
