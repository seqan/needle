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
outputs:
  []
label: needle-minimiser
doc: ""
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - needle
  - minimiser
