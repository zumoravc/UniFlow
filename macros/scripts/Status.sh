#!/bin/bash

name="vpacik"

# vojtech.pacik@cern.ch

tmp1="jobs-all.tmp"
tmp2="jobs-sel.tmp"
cat /dev/null > $tmp1
cat /dev/null > $tmp2

alien_top -user $name -all_status > $tmp1

cat $tmp1 | grep "pcapiserv" > $tmp2; n1=$(wc -l $tmp2 | awk 'NR==1{print $1;}')
cat $tmp1 | grep "pcapiserv" | grep "DONE" > $tmp2; n2=$(wc -l $tmp2 | awk 'NR==1{print $1;}')
echo "-------------------------------"
echo "Master jobs:            "$n1
echo "Done master jobs:       "$n2
echo "Master jobs completed:  "$(echo -e "if $n1 > 0: print '%.1f' %(100.*$n2/$n1) \nelse: print 0" | python 2>/dev/null)" %"
echo "-------------------------------"

cat $tmp1 | grep "aliendb" > $tmp2; n1=$(wc -l $tmp2 | awk 'NR==1{print $1;}')
cat $tmp1 | grep "aliendb" | grep "ERROR" > $tmp2
cat $tmp1 | grep "aliendb" | grep "ZOMBIE" >> $tmp2
cat $tmp1 | grep "aliendb" | grep "EXPIRED" >> $tmp2
n3=$(wc -l $tmp2 | awk 'NR==1{print $1;}')
cat $tmp1 | grep "aliendb" | grep "DONE" > $tmp2; n2=$(wc -l $tmp2 | awk 'NR==1{print $1;}')

echo "Jobs:             "$n1
echo "Errors:           "$n3
echo "Unfinished jobs:  "$(expr $n1 - $n2 - $n3)
echo "Done jobs:        "$n2
echo "Jobs completed:   "$(echo -e "if $n1 > 0: print '%.1f' %(100.*$n2/$n1) \nelse: print 0" | python 2>/dev/null)" %"
echo "-------------------------------"
date

wait
rm -f $tmp1 2>/dev/null
rm -f $tmp2 2>/dev/null
