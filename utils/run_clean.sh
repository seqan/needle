#!/bin/bash
set -eu

metagraph="Set path to metagraph executable"
# Note: You need to have a dbg_files.txt, which is a file where in each line there is a path to the dbg files with the ending ".dbg" (created with run_build.sh)

mkdir single_dbgs/clean/

while read -r filename; do
	$metagraph clean -p4 --to-fasta --primary-kmers --smoothing-window 1 -o single_dbgs/clean/$(basename $filename)  --count-kmers --count-width 32  $filename
done  < dbg_files.txt
