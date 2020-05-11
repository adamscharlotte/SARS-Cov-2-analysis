cwlVersion: v1.0
class: CommandLineTool

baseCommand: [Rscript]

requirements:
 - class: InlineJavascriptRequirement

inputs:
    pia_genes_R:
        type: File
        inputBinding:
            position: 1

    proteins_csv:
        type: File
        inputBinding:
            position: 2

    bait:
        type: string
        inputBinding:
            position: 3

outputs:
    genes_csv:
        type: File
        outputBinding: 
            glob: $(inputs.bait).csv

#    ANN-SoLo_identifications:
 #       type: stdout
#stdout: tibble.txt          #this will appear in the working directory
