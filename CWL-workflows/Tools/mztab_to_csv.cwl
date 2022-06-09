cwlVersion: v1.0
class: CommandLineTool

baseCommand: [python]

requirements:
 - class: InlineJavascriptRequirement

inputs:
    mztab_to_csv_py:
        type: File
        inputBinding:
            position: 1

    mztab:
        type: File
        inputBinding:
            position: 2

    filename:
        type: string
        inputBinding:
            position: 3

outputs:
    csv:
        type: File
        outputBinding: 
            glob: $(inputs.filename)*.csv
#    tibbles:
 #       type: stdout
#stdout: tibble.txt          #this will appear in the working directory
