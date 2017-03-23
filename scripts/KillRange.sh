#!/bin/bash

name="vpacik"


# jaroslav.adam@cern.ch

tmpfile="myjobs.tmp"
cat /dev/null > $tmpfile

alien_top -user $name -all_status | grep "pcapiserv" > $tmpfile

while read i
do
  id=`echo $i | awk 'NR==1{print $1;}'`
  #echo $id

  if [[ $id -gt 779532835 ]]; then
  	alien_kill $id
  fi


done < $tmpfile

rm -f $tmpfile 2>/dev/null

