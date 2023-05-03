inputs:
  - doc: "Define k-mer size for the minimisers. Default: 20. "
    id: kmer
    type:
      - "null"
      - long
    inputBinding:
      prefix: --kmer
  - doc: "Define window size for the minimisers. Default: 60. Default: 0. "
    id: window
    type:
      - "null"
      - long
    inputBinding:
      prefix: --window
  - doc: "Define a shape for the minimisers by the decimal of a bitvector, where 0 symbolizes a position to be ignored, 1 a position considered. Default: ungapped. Default: 0. "
    id: shape
    type:
      - "null"
      - long
    inputBinding:
      prefix: --shape
  - doc: "Define seed for the minimisers. Default: 0. "
    id: seed
    type:
      - "null"
      - long
    inputBinding:
      prefix: --seed
  - doc: "Directory, where output files should be saved. Default: \"./\". "
    id: out
    type:
      - "null"
      - string
    inputBinding:
      prefix: --out
  - doc: "Number of threads to use. Default: 1. "
    id: threads
    type:
      - "null"
      - long
    inputBinding:
      prefix: --threads
  - doc: "List of bin false positive rate per expression level. If only one is given, then that fpr is used for all expression levels. Default: []. "
    id: fpr
    type:
      - "null"
      - boolean
    inputBinding:
      prefix: --fpr
  - doc: "Which expression thresholds should be used for constructing the IBFs. Default: []. "
    id: expression_thresholds
    type:
      - "null"
      - boolean
    inputBinding:
      prefix: --expression_thresholds
  - doc: "Number of expression thresholds. Can be set alternatively to expression_thresholds, then the expression thresholds are determined automatically. Default: 0. "
    id: number_expression_thresholds
    type:
      - "null"
      - long
    inputBinding:
      prefix: --number_expression_thresholds
  - doc: "Number of hash functions that should be used when constructing one IBF. Default: 1. "
    id: hash
    type:
      - "null"
      - long
    inputBinding:
      prefix: --hash
  - doc: "Sequence file containing minimizers, only those minimizers will be considered. Default: \"\". "
    id: include
    type:
      - "null"
      - string
    inputBinding:
      prefix: --include
  - doc: "Sequence file containing minimizers that should not be stored. Default: \"\". "
    id: exclude
    type:
      - "null"
      - string
    inputBinding:
      prefix: --exclude
  - doc: "Define which samples belong together, sum has to be equal to number of sequence files. Default: Every sequence file is one sample from one experiment. Default: []. "
    id: samples
    type:
      - "null"
      - boolean
    inputBinding:
      prefix: --samples
  - doc: "Define for each sample, what number of found minimisers should be considered the result of a sequencing error and therefore be ignored. Default: Every sample has an automatically generated cutoff, which is based on the file size. Default: []. "
    id: cutoff
    type:
      - "null"
      - boolean
    inputBinding:
      prefix: --cutoff
  - doc: "If set, names of the experiments are stored in a txt file. Default: 0. "
    id: experiment-names
    type:
      - "null"
      - boolean
    inputBinding:
      prefix: --experiment-names
  - doc: "Sequence file containing minimizers, only those minimizers will be considered for determining the expression thresholds. Default: \"\". "
    id: levels-by-genome
    type:
      - "null"
      - string
    inputBinding:
      prefix: --levels-by-genome
outputs:
  []
label: needle-ibf
doc: ""
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - needle
  - ibf
