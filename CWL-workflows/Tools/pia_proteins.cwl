cwlVersion: v1.0
class: CommandLineTool

baseCommand: [Rscript]

requirements:
 - class: InlineJavascriptRequirement

inputs:
    pia_proteins_R:
        type: File
        inputBinding:
            position: 1

    csv:
        type: File
        inputBinding:
            position: 2

    filename:
        type: string
        inputBinding:
            position: 3

outputs:
    proteins_csv:
        type: File
        outputBinding: 
            glob: $(inputs.filename).csv

#    ANN-SoLo_identifications:
 #       type: stdout
#stdout: tibble.txt          #this will appear in the working directory
