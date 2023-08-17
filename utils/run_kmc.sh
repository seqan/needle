#!/bin/bash
set -eu

kmc="Path to kmc executable"

mkdir kmc_files
mkdir kmc_tmp
while read -r filename threshold; do
	$kmc -t64 -r -k21 -ci$threshold -cs65535 -hp $filename kmc_files/$(basename $filename) kmc_tmp/
	echo $filename
done  < samples.in
