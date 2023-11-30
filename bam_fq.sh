#!/usr/bin/bash

currpath=$(pwd)
ls -l *.bam.1 |awk '{print $9}' | while read i
do 
bamtofastq --nthreads=32 \
  $currpath/$i \
$currpath/fastq_from_bam/$i
done
