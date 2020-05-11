cwlVersion: v1.0
class: CommandLineTool

baseCommand: [python3]

requirements:
 - class: InlineJavascriptRequirement

inputs:
    mztab_to_idxml_py:
        type: File
        inputBinding:
            position: 1

    mztab:
        type: File
        inputBinding:
            position: 2

    bait:
        type: string
        inputBinding:
            position: 3

outputs:
    idxml:
        type: File
        outputBinding: 
            glob: $(inputs.bait).idxml
#    tibbles:
 #       type: stdout
#stdout: tibble.txt          #this will appear in the working directory
