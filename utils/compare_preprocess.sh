#!/bin/bash

needle="Set path to needle executable"
input_dir="Set to directory that contains the sequence files"

# Comparison for 4 Threads, change thread number in each step to obtain numbers for 4 threads
/usr/bin/time -v -o needle_21_21_preprocess4.time $needle minimiser -k 21 -w 21 -t 4 --cutoff 49 $input_dir/SRR1313229.fastq.gz $input_dir/SRR1313228.fastq.gz  $input_dir/SRR1313227.fastq.gz $input_dir/SRR1313226.fastq.gz
/usr/bin/time -v -o needle_25_21_preprocess4.time $needle minimiser -k 21 -w 25 -t 4 --cutoff 49 $input_dir/SRR1313229.fastq.gz $input_dir/SRR1313228.fastq.gz  $input_dir/SRR1313227.fastq.gz $input_dir/SRR1313226.fastq.gz
/usr/bin/time -v -o needle_21_41_preprocess4.time $needle minimiser -k 21 -w 41 -t 4 --cutoff 49 $input_dir/SRR1313229.fastq.gz $input_dir/SRR1313228.fastq.gz  $input_dir/SRR1313227.fastq.gz $input_dir/SRR1313226.fastq.gz

# For both bcalm and kmc a file named "files.lst" with the path to the four files needs to be created.
/usr/bin/time -v -o bcalm_preprocess4.time bash run_bcalm4.sh
/usr/bin/time -v -o kmc_preprocess4.time bash run_kmc4.sh
