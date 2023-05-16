import sys

import numpy as np

infile = sys.argv[1]
expfile = sys.argv[2]
outfile = sys.argv[3]

genes = {}
with open(infile, 'r') as f:
    for line in f:
        gene = line.split('|')[5]
        expressions = [int(x) for x in line.strip().split()[1:]]
        if (gene in genes):
            genes[gene].append(expressions)
        else:
            genes.update({gene:[expressions]})

expnames = []
with open(expfile, 'r') as f:
    for line in f:
        expnames.append(line.strip())

with open(outfile, 'w') as o:
    o.write("gene\t")
    for i in range(len(expnames) - 1):
        o.write(expnames[i])
        o.write("\t")
    o.write(expnames[len(expnames) - 1])
    o.write("\n")
    for gene in genes.keys():
        o.write(gene)
        o.write("\t")
        mean_exp = np.array(genes[gene])
        #print(mean_exp)
        mean_exp = np.mean(mean_exp, axis = 0)
        #print(gene, mean_exp)
        for k in range(len(mean_exp) - 1):
            o.write(str(int(mean_exp[k])))
            o.write("\t")
        o.write(str(int(mean_exp[len(mean_exp) - 1])))
        o.write("\n")
