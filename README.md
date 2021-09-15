## Needle
Needle provides a space-efficient data structure to index a large amount of NGS data and allows fast searches through these indices.
Due to the space-efficiency of one index, it is affordable to create multiple indices with different expression rates. Therefore, a semi-quantitative analysis of the data becomes possible. Needle is based on Interleaved Bloom Filters, which is a compact and efficient structure to store multiple Bloom Filters. Furthermore, Needle uses a windowing scheme (also called minimisers) to reduce the amount of data to store.

## Build

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

## Create an IBF
In order to create an IBF a number of sequence files have to be given. All sequence file formats from seqan3 are accepted as an input (fasta, fastq, embl,... and their compressed forms). With the parameter m can be defined, which of these sequence files belong together, either because they are the result of paired-end sequencing or they are multiple replicates of the same experiment. If no specification with m is given, every sequence file is seen as one experiment. For paired-end experiments one can use the flag '--paired' to indicate this, so two consecutive sequence files are seen as belonging together. (This is equivalent to using -m 2 for all experiments.)
Besides, the false positive rate of the IBF has to be specified with parameter f.
Use -h/--help for more information and to see further parameters.

The following example creates an IBF for two experiments for the expression levels 4 and 32. Both experiments had two replicates, therefore m is used to specify this. With c a compressed IBF is created.

```
./bin/needle ibf ../needle/test/data/exp_*.fasta --samples 2 --samples 2 -e 4 -e 32 -f 0.3 -c -o example

// Or with flag paired
./bin/needle ibf ../needle/test/data/exp_*.fasta --paired -e 16 -e 32 -f 0.3 -c -o example
```

## Calculate Minimisers
In case one is only interested in the minimisers or wants to preprocess the data first before creating an IBF, the function minimiser can be used. It calculates the minimisers of given experiments and stores their hash values and their occurrences in a binary file named ".minimiser". Furthermore, a txt file is created where all used arguments are stored (like used k-mer size or window size), the used expression levels and the minimiser counts per expression level.

The following command calculates the minimisers in the two experiments.
```
./bin/needle minimiser ../needle/test/data/exp_*.fasta -samples 2 -samples 2
```

A minimiser file is a binary file containing the following data:
- number of minimisers (uint64_t)
- kmer-size (uint8_t)
- window-size (uint32_t)
- seed (uint64_t)
- flag which is true, if shape is ungapped (bool)
- shape (uint64_t), if flag is false
- all minimiser hashes (uint64_t) with their occurrences (uint16_t)

Based on a minimiser file the ibfs can be computed by using the following command:
```
./bin/needle ibfmin exp*.minimiser -e 16 -e 32  -f 0.3 -c -o example
```

## Estimate
To estimate the expression value of one transcript a sequence file has to be given. Use the parameter "-i" to define where the IBFs can be found (should be equal with "-o" in the previous commands).
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
