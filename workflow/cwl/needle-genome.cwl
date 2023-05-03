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
  - doc: "Please provide one sequence file with minimizers to ignore. Default: \"\". "
    id: exclude
    type:
      - "null"
      - string
    inputBinding:
      prefix: --exclude
outputs:
  []
label: needle-genome
doc: ""
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - needle
  - genome
