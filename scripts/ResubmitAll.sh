#!/bin/bash
name="vpacik"

tmpfil1="myjobs.tmp"
tmpfil2="errjobs.tmp"
#cat /dev/null > $tmpfil1
#cat /dev/null > $tmpfil2

alien_top -user $name -all_status > $tmpfil1

cat $tmpfil1 | grep "ERROR" >> $tmpfil2
cat $tmpfil1 | grep "ZOMBIE" >> $tmpfil2
cat $tmpfil1 | grep "EXPIRED" >> $tmpfil2

while read i
do
  id=`echo $i | awk 'NR==1{print $1;}'`
  #echo $id

  alien_resubmit $id

done < $tmpfil2

#wait
#rm -f $tmpfil1 2>/dev/null
# rm -f $tmpfil2 2>/dev/null

