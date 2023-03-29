# Use: python3 unify_results.py [output_file] [input_files in one file] [file with accession associated with condition]
# python3 unify_results.py deseq_in.csv infiles.lst accession.lst
import os
import sys

outfile = sys.argv[1]
infiles = sys.argv[2]
cancer = sys.argv[3]

def get_gene_exp(inf, outf):
    genes = {}
    with open(inf, 'r') as f:
        for line in f:
            gene = line.strip().split('|')[5]
            val = int(line.strip().split()[1])
            pc = line.strip().split('|')[7]
            if (pc == "protein_coding"):
                if (gene in genes):
                    genes[gene].append(val)
                else:
                    genes.update({gene: [val]})

    with open(outf, 'w') as o:
        for gene in sorted(genes.keys()):
            o.write(gene)
            o.write('\t')
            o.write(str(round(sum(genes[gene])*1.0/len(genes[gene]))))
            o.write('\n')

files = []

with open(infiles, 'r') as f:
    for line in f:
        files.append(line.strip())

conditions = []

with open(cancer, 'r') as f:
    for line in f:
        if (line.strip().split()[0] != "condition"):
            conditions.append(line.strip().split()[0])

counts = []
genes = []
i = 0
for file in files:
    counts.append([])
    if not (os.path.isfile("Exp_"+file)):
        get_gene_exp(file, "Exp_"+file)
    with open("Exp_"+file, 'r') as f:
        for line in f:
            counts[i].append(int(line.strip().split()[1]))
            if (i == 0):
                genes.append(line.strip().split()[0])
    i += 1

with open(outfile, 'w') as o:
    o.write("gene id\t")
    for j in range(i):
        o.write(conditions[j])
        if (j < len(counts) - 1):
            o.write('\t')
    o.write('\n')
    for k in range(len(genes)):
        o.write(genes[k])
        o.write('\t')

        for j in range(i):
            if (k == 0):
                print(j, counts[j][k], (j < len(counts) - 1))
            o.write(str(counts[j][k]))
            if (j < len(counts) - 1):
                o.write('\t')
            elif (k == 0):
                print(j)
        o.write('\n')
