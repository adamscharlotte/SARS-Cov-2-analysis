cwlVersion: v1.0
class: CommandLineTool

baseCommand: [Rscript]

requirements:
 - class: InlineJavascriptRequirement
inputs:
    modifications_R:
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

    mass_tolerance:
        type: string
        inputBinding:
            position: 4

outputs:
    tol_py_csv:
        type: File
        outputBinding: 
            glob: $(inputs.bait)_tol_py.csv
#    tibbles:
 #       type: stdout
#stdout: tibble.txt          #this will appear in the working directory