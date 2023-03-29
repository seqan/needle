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
