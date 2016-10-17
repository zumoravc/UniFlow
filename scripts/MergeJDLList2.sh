#!/bin/bash

runList="10eAOD135test.txt"
prefix="000"

echo "### Merging jobs via JDL from $runList ###"
for i in `cat $runList`
do
  echo "=== Merging $prefix$i ==="
  path="test/$prefix$i"
  cd $path/
  root -l -b -q "~/Codes/ALICE/StrangenessInJets/macros/runEMCalJetGrid.C(kTRUE,kTRUE,\"terminate\",\"$i\")" > stdMerge.log 2>&1 &
  cd ../../
  sleep 20
done
echo "### Merging via JDL DONE ###"
