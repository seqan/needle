## Needle
Needle provides a space-efficient data structure to index a large amount of NGS data and allows fast searches through these indices.
Due to the space-efficiency of one index, it is affordable to create multiple indices with different expression rates. Therefore, a semi-quantitative analysis of the data becomes possible. Needle is based on Interleaved Bloom Filters, which is a compact and efficient structure to store multiple Bloom Filters. Furthermore, Needle uses a windowing scheme (also called Minimizers) to reduce the amount of data to store.  

## Build

Needle can be built by following these commands:

```
git clone --recurse-submodules https://github.com/MitraDarja/needle.git
mkdir build-needle && cd build-needle
cmake ../needle
make
```

Run test to check, if Needle is working as intended. All tests should pass.

```
./test/api/test-needle
```

## Create an IBF
In order to create an IBF a number of sequence files have to be given. All sequence file formats from seqan3 are accepted as an input (fasta, fastq, embl,... and their compressed forms). With the parameter m can be defined, which of these sequence files belong together, either because they are the result of paired-end sequencing or they are multiple replicates of the same experiment. If no specification with m is given, every sequence file is seen as one experiment. For paired-end experiments one can use the flag q to indicate this, so two consecutive sequence files are seen as belonging together. (This is equivalent to using -m 2 for all experiments.)
Besides, the sequence file the bin size of the IBF has to be specified with parameter b. Good sizes for Bloom Filters for one experiment can be calculated with this [calculator] (https://hur.st/bloomfilter/?n=&p=5.0E-2&m=6559922&k=1).
Use -h/--help for more information and to see further parameters.

The following example creates an IBF for two experiments for the expression rate 0.5. Both experiments had two replicates, therefore m is used to specify this. With c a compressed IBF is created.

```
./needle ibf ../needle/test/data/exp_*.fasta -m 2 -m 2 -e 0.5 -b 100000 -c

// Or with flag q
./needle ibf ../needle/test/data/exp_*.fasta -q -e 0.5 -b 100000 -c
```

## Calculate Minimizers
In case one is only interested in the minimizers or wants to preprocess the data first before creating an IBF, the function minimizer can be used. It calculates the minimizers of given experiments and stores their hash values and their occurences in a binary file named ".minimizer". Furthermore, a txt file is created where all used arguments are stored (like used k-mer size or window size), the used expression levels and the minimizer counts per expression level.

The following command calculates the minimizers in the two experiments 0 and 1 for three different expression levels.
```
./needle minimizer ../needle/test/data/exp_*.fasta -m 2 -m 2 -e 0 -e 1 -e 4
```

A header file experiment name ".header" looks like this:
```
10322096095657499358 20 60 0 median 29 // seed k-mer_size window_size shape normalization_method normalized_expression_value
0 1 4                                  // expression levels
62496 6116 25                          // minimizer count per expression level, so 62496 for 0, 6116 for 1, 25 for 4
```

With the function stats some information about the counts on different expression levels and the normalized expression values can be calculated. (As shown in the following example.)

```
./needle stats Header_exp_01.txt Header_exp_11.txt

```

## Search
For a search at least one transcript to be searched in a sequence file format has to be given and an expression rate. If a compressed IBF was created the search needs the information that it is a compressed IBF as well.
Use -h/--help for more information and to see further parameters.
The following example searches for one gene, which is expressed in the first experiment at expression rate 0.25 and in the second at expression rate 1. Therefore, it should be found only in the second experiment but not the first.

```
./needle search ../needle/test/data/gene.fasta -e 0.5 -c
```
This results in:
```
IBF_0.500000
Results: 0 1
```

## Insert
After an IBF is created, it is possible to add further experiments with the insert function, but only if it is an uncompressed IBF. In order to have a consistent IBF the same parameters (k-mer size, window size, normalization method, ...) should be chosen.

```
./needle ibf ../needle/test/data/exp_0*.fasta -m 2 -e 0.5 -b 100000 // Create IBF with one experiment
./needle insert ../needle/test/data/exp_1*.fasta -m 2 -e 0.5        // Adds the second experiment
```
