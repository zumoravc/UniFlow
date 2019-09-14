#!/bin/bash

list="killList.txt"

for i in `cat $list`
do
  alien_kill $i
done
