#!/bin/bash

### Pro vycisteni spatneho / neplanovaneho mergeViaJDL


runList="10eAOD135test.txt"
prefix="000"

for i in `cat $runList`
do
  echo "=== $i ==="
  alien_rm /alice/cern.ch/user/v/vpacik/work/EMCalV0sAnomaly10e_run2_$i/output/000/Stage_*.xml
  alien_rm /alice/cern.ch/user/v/vpacik/work/EMCalV0sAnomaly10e_run2_$i/output/000/*.root
  alien_rm /alice/cern.ch/user/v/vpacik/work/EMCalV0sAnomaly10e_run2_$i/output/000/*.zip
  alien_rmdir /alice/cern.ch/user/v/vpacik/work/EMCalV0sAnomaly10e_run2_$i/output/000/Stage_*
	
done
