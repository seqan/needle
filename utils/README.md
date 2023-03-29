### Evaluation of Needle count

# Differential Expression
Download the sequencing experiments listed in true_breast.csv and the human transcripts from gencode as a fasta file.
Then create a Needle-genome file via:

```
needle genome -k 19 -w 19 -o 19_19_ gencode.v43.transcripts.fa
```

Afterwards, Needle count is run via:

```
needle count --genome 19_19_gencode.v43.transcripts.genome -k 19 -w 19 -o 19_19_ --include gencode.v43.transcripts.fa SRR66711*.gz
```

To combine all results in one file, use the provided python script.

```
python3 unify_results2.py deseq_in.csv infiles.lst acc_info.lst
```

Now, you can run the RScript "" to analyze the differential expressed genes and generate the heatmap for the 67 
differential expressed genes.
