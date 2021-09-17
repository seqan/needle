## Needle

[![Build Status](https://github.com/seqan/app-template/workflows/App%20CI/badge.svg)](https://github.com/seqan/needle/actions?query=branch%3Amaster+workflow%3A%22App+CI%22)

Needle is a tool for semi-quantitative analysis of very large collections of nucleotide sequences.
Needle stores its data in multiple interleaved Bloom filter, a fast and space efficient probabilistic data structure and uses a windowing scheme (also called minimisers) to reduce the amount of data to store. How many interleaved Bloom filter are used is defined by the user. Each interleaved Bloom filter has a so called expression threshold and stores minimisers with an occurrence greater than or equal to its own expression threshold and smaller than the next biggest expression threshold (if there is no bigger expression threshold, all greater than or equal to the threshold are stored). These expression thresholds are then used during the query (called estimate) to approximate the expression values of given transcripts.

## Install

Needle can be built by following these commands:

```
git clone --recurse-submodules https://github.com/seqan/needle.git
mkdir build-needle && cd build-needle
cmake ../needle
make
```

Run test to check, if Needle is working as intended. All tests should pass.

```
make test
```

If you are interested in building the documentation, just use the command: `make doc`

## Build
In order to build a Needle index a number of sequence files have to be given. All sequence file formats supported by seqan3 are accepted as an input (fasta, fastq, embl,... and their compressed forms). The flag `--paired` in the example below indicates that the given sequence files are paired-end experiments. Furthermore, the false positive rate has to be specified with the parameter `f`.
Use -h/--help for more information and to see further parameters. The flag `-c` can be used to build a compressed Needle index.

The following example creates a compressed Needle index for two paired-end experiments for the expression thresholds 4 and 32.

```
./bin/needle ibf ../needle/test/data/exp_*.fasta --paired -e 16 -e 32 -f 0.3 -c -o example
```

Although, this works. It is recommended to calculate the minimisers beforehand by using the option `minimisers`. It calculates the minimisers of given experiments and stores their hash values and their occurrences in a binary file named ".minimiser".

The following command calculates the minimisers in the two experiments.
```
./bin/needle minimiser ../needle/test/data/exp_*.fasta --paired
```

A minimiser file is a binary file containing the following data:
- number of minimisers (uint64_t)
- kmer-size (uint8_t)
- window-size (uint32_t)
- seed (uint64_t)
- flag which is true, if shape is ungapped (bool)
- shape (uint64_t), if flag is false
- all minimiser hashes (uint64_t) with their occurrences (uint16_t)

Based on the minimiser files the Needle index can be computed by using the following command:
```
./bin/needle ibfmin exp*.minimiser -e 16 -e 32  -f 0.3 -c -o example
```

## Estimate
To estimate the expression value of one transcript a sequence file has to be given. Use the parameter "-i" to define where the Needle index can be found (should be equal with "-o" in the previous commands).
Use -h/--help for more information and to see further parameters.
The following example searches for one gene, which is expressed in the first experiment with expression 6 and in the second with expression 37. Therefore, it should be found only in the second experiment but not the first when using expression levels of 16 and 32.

```
./bin/needle estimate ../needle/test/data/gene.fasta -i example
```

The created file "expressions.out" (if you prefer a different name, use "-o") should contain the following:
```
GeneA   0      32
```

## Note

This app was created with the [seqan3 app-template](https://github.com/seqan/app-template).
