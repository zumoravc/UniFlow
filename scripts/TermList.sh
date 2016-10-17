#!/bin/bash

runList="10eAOD135full.txt"
prefix="000"

echo "### Terminating & downloading jobs from $runList ###"
for i in `cat $runList`
do
  echo "=== Terminating $prefix$i ==="
  path="LHC10e/$prefix$i"
  cd $path/
  root -l -b -q "runEMCalJetGrid.C(kTRUE,kFALSE,\"terminate\",\"$i\")" > stdTerm.log 2>&1 &
  cd ../../
  sleep 20
done
echo "### Terminating & downloading DONE ###"
