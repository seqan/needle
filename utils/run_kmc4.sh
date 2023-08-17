#!/bin/bash

kmc="Set path to kmc executable"

#get fastq.gz and launch bcalm on each file
while read -r filename; do
    $kmc -t4 -r -k21 -ci50 -cs65535 -hp -m500 $filename  $(basename $filename .res) kmc_tmp/
done < files.lst
