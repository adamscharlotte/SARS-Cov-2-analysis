cwlVersion: v1.0
class: CommandLineTool

baseCommand: [python]

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

    idxml:
        type: string
        inputBinding:
            position: 3

outputs:
    idxml:
        type: File
        outputBinding: 
            glob: $(inputs.idxml)
#    tibbles:
 #       type: stdout
#stdout: tibble.txt          #this will appear in the working directory
