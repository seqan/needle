# Evaluation of Needle count

## Differential Expression
Download the sequencing experiments listed in accession.lst and the human transcripts from gencode as a fasta file.
Then create a Needle-genome file via:

```
needle genome -k 19 -w 19 -o 19_19_ gencode.fa
needle genome -k 19 -w 23 -o 19_23_ gencode.fa
needle genome -k 19 -w 39 -o 19_39_ gencode.fa
```

Afterwards, Needle count is run via:

```
needle count --genome 19_19_gencode.genome -k 19 -w 19 -o 19_19_ --include gencode.fa SRR66711*.gz
needle count --genome 19_23_gencode.genome -k 19 -w 23 -o 19_23_ --include gencode.fa SRR66711*.gz
needle count --genome 19_39_gencode.genome -k 19 -w 39 -o 19_39_ --include gencode.fa SRR66711*.gz
```

To combine all results in one file, use the provided python script. The infiles.lst can be created by storing the output files of the count command.

```
ls 19_19_*count.out > 19_19_infiles.lst
ls 19_23_*count.out > 19_23_infiles.lst
ls 19_39_*count.out > 19_39_infiles.lst
python3 unify_results.py deseq_in.csv 19_19_infiles.lst accession.lst
python3 unify_results.py deseq_in_23.csv 19_23_infiles.lst accession.lst
python3 unify_results.py deseq_in_39.csv 19_39_infiles.lst accession.lst
```

Now, you can run the RScript "breastcancer.R" to analyze the differential expressed genes and generate the heatmap for the 67 
differential expressed genes from https://doi.org/10.1016/j.dib.2018.03.079.


# Evaluation of Needle

## Parallelization

In order to compare the two parallelization methods of Needle when computing minimizers, four sequencing experiments were tested.
For the simple version the following command was used with either `X=1` or `X=4` threads.
```
needle minimiser -k 21 -w 25 -t X SRR1313229.fastq.gz SRR1313228.fastq.gz  SRR1313227.fastq.gz SRR1313226.fastq.gz
```
For the ram friendly version a `--ram` was added to the previous command.

## Differential Expression
Download the sequencing experiments of the GEO experiment with the accession number GSE58135 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58135).
Then create a Needle index in the following way, assuming all sequencing files are stored in a folder named GSE58135:

```
# Create minimiser files
minimiser -k 20 -w 24 -o Cutoff_1_ --cutoff 1 GSE58135/*.fastq.gz -t 4
# Create index
needle ibfmin -o 24_20_Cutoff_1_Compressed_ -l 15 $(ls -v Cutoff_1_*minimiser) -f 0.05 -c
```

Afterwards, download a gencode file to determine the quantification of gene transcripts and use the provided python script afterwards to get the expression per gene.

```
needle estimate -i 24_20_Cutoff_1_Compressed_ gencode.fa.gz -o Cutoff_1_Compressed_expressions_all.out
python3 needle_gene_exp.py Cutoff_1_Compressed_expressions_all.out experiments.out Cutoff_1_Compressed_20_24_gene_expression_all.csv
```

Now, you can run the RScript "breastcancer_big.R" to analye the differential expressed genes for different breast cancer gene signatures.

## Gene Ontology Analysis

The gene ontology analysis can be repeated by querying all proteincoding human transcripts in the Needle index created [here](https://github.com/MitraDarja/analysis_needle/blob/main/run_large_dataset.sh) by running:
```
needle estimate -i w_21/SRR_Compressed gencode_all_human_pc.fa.gz -o gencode_srr.out
```

Then obtain all expressions per gene by using the file [sras_1742.lst](https://github.com/MitraDarja/analysis_needle/blob/main/data/sras_1742.lst):
```
python3 needle_gene_exp.py gencode_srr.out sras_1742.lst gencode_srr.csv
```

Afterwards, the differential expressed genes per tissue type can be obtained by:
```
RScrip tissue_type.R
```

The differential expressed genes can then by obtained via ShinyGo using Ontology.Jensen.TISSUES as a database.

