## Needle
Needle provides a space-efficient data structure to index a large amount of NGS data and allows fast searches through these indices.
Due to the space-efficiency of one index, it is affordable to create multiple indices with different expression rates. Therefore, a semi-quantitative analysis of the data becomes possible. Needle is based on Interleaved Bloom Filters, which is a compact and efficient structure to store multiple Bloom Filters. Furthermore, Needle uses a windowing scheme (also called Minimizers) to reduce the amount of data to store.  

## Build

Needle is depending on the seqan3 library (https://github.com/seqan/seqan3), at the moment it is necessary to use Enrico Seiler's branch "feature/binning_directory" (https://github.com/eseiler/seqan3/tree/feature/binning_directory), where the IBF is implemented. Soon, this branch should be included in the seqan3 library.
Assuming seqan3 can be found in "${CMAKE_SOURCE_DIR}/../", Needle can be build following these commands:

```
git clone https://github.com/MitraDarja/needle.git
mkdir build-needle && cd build-needle
cmake ../needle
make
```

## Create an IBF
In order to create an IBF a number of sequence files have to be given. All sequence file formats from seqan3 are accepted as an input (fasta, fastq, embl,... and their compressed forms). With the parameter m can be defined, which of these sequence files belong together, either because they are the result of paired-end sequencing or they are multiple replicates of the same experiment. If no specification with m is given, every sequence file is seen as one experiment.
Besides, the sequence file a size of the IBF has to be specified with parameter l. Good sizes for Bloom Filters for one experiment can be calculated with this calculator (https://hur.st/bloomfilter/?n=&p=5.0E-2&m=6559922&k=1), which should be then multiplied with the number of given experiments.
Use -h/--help for more information and to see further parameters.

The following example creates an IBF for two experiments for the expression rate 0.5. Both experiments had two replicates, therefore m is used to specify this. With c a compressed IBF is created.

```
./needle-ibf ../../needle/example/exp_*.fasta -m 2 -m 2 -e 0.5 -l 1559922 -c
```

## Search
For a search at least one transcript to be searched in a sequence file format has to be given and an expression rate. If a compressed IBF was created the search needs the information that it is a compressed IBF as well.
Use -h/--help for more information and to see further parameters.
The following example searches for one gene, which is expressed in the first experiment at expression rate 0.25 and in the second at expression rate 1. Therefore, it should be found only in the second experiment but not the first.

```
./needle-search ../../needle/example/gene.fasta -e 0.5 -c
```
This results in:
```
IBF_0.500000
Results: [0,1]
```
