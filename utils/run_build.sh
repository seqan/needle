#!/bin/bash
set -eu

metagraph="Set path to metagraph executable"
# Note: You need to have a kmc_files.lst, which is a file where in each line there is a path to the kmc_files with the ending ".kmc_suf"

mkdir single_dbgs/
while read -r filename; do
	$metagraph build --state fast --mode canonical --parallel 16 --count-kmers --count-width 32 -k 21 --mem-cap-gb 8 -o single_dbgs/$(basename ${filename})   $filename
	echo $filename
done  < kmc_files.lst
