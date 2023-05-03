inputs:
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
  - doc: "Sequence file containing minimizers, only those minimizers will be considered for determining the expression thresholds. Default: \"\". "
    id: levels-by-genome
    type:
      - "null"
      - string
    inputBinding:
      prefix: --levels-by-genome
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
outputs:
  []
label: needle-ibfmin
doc: ""
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - needle
  - ibfmin
