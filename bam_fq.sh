#!/usr/bin/bash


#reference to https://github.com/10XGenomics/bamtofastq
currpath=$(pwd)
ls -l *.bam.1 |awk '{print $9}' | while read i
do 
bamtofastq --nthreads=32 \
  $currpath/$i \
$currpath/fastq_from_bam/$i
done
