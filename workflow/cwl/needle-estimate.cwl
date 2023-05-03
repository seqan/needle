inputs:
  - doc: "Directory where input files can be found. Default: \"./\". "
    id: in
    type:
      - "null"
      - string
    inputBinding:
      prefix: --in
  - doc: "Directory, where output files should be saved. Default: \"expressions.out\". "
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
outputs:
  []
label: needle-estimate
doc: ""
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - needle
  - estimate
