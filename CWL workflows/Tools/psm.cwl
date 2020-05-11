cwlVersion: v1.0
class: CommandLineTool

baseCommand: [Rscript]

requirements:
 - class: InlineJavascriptRequirement

inputs:
    psm_R:
        type: File
        inputBinding:
            position: 1

    csv:
        type: File
        inputBinding:
            position: 2

    bait:
        type: string
        inputBinding:
            position: 3

outputs:
    csv:
        type: File
        outputBinding: 
            glob: $(inputs.bait).csv