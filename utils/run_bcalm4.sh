#!/bin/bash
# Based on https://github.com/kamimrcht/REINDEER/blob/master/reproduce_manuscript_results/bcalm_2585.sh

bcalm="Set path to bcalm2 executable"

#get fastq.gz and launch bcalm on each file
while read -r filename; do
    $bcalm -in $filename -kmer-size 21 -abundance-min 50 -nb-cores 2 -max-memory 500000 -out-tmp bcalm2 -out-dir bcalm2
done < files.lst
