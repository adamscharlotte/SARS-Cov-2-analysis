cwlVersion: v1.0
class: CommandLineTool

baseCommand: [python3]

requirements:
 - class: InlineJavascriptRequirement

inputs:
    modifications_faster_py:
        type: File
        inputBinding:
            position: 1

    bait:
        type: string
        inputBinding:
            position: 2

    unimod_path:
        type: string
        inputBinding:
            position: 3
 
    tol_py_csv:
        type: File
        inputBinding:
            position: 4

outputs:
    mod_csv:
        type: File
        outputBinding: 
            glob: $(inputs.bait).csv
#    tibbles:
 #       type: stdout
#stdout: tibble.txt          #this will appear in the working directory