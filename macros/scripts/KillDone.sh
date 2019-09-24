#!/bin/bash

name="zumoravc"


# zuzana.moravcova@cern.ch

tmpfile="myjobs.tmp"
cat /dev/null > $tmpfile

alien_top -user $name -all_status | grep "aliendb" | grep "DONE" > $tmpfile

while read i
do
  id=`echo $i | awk 'NR==1{print $1;}'`
  #echo $id

  alien_kill $id

done < $tmpfile

rm -f $tmpfile 2>/dev/null
