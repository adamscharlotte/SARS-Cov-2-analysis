cwlVersion: v1.0
class: CommandLineTool

baseCommand: [java]

requirements:
 - class: InlineJavascriptRequirement

inputs:
    cp:
        type: File
        inputBinding:
            position: 1
            prefix: -cp

    compiler:
        type: string
        inputBinding:
            position: 2

    idxml:
        type: File
        inputBinding:
            position: 3
            prefix: -infile

    bait:
        type: string
        inputBinding:
            position: 4    

    xml_output:
        type: string
        inputBinding:
            position: 5
            prefix: -outfile

    name:
        type: string
        inputBinding:
            position: 6
            prefix: -name

outputs:
    xml:
        type: File
        outputBinding: 
            glob: $(inputs.bait).xml
#    tibbles:
 #       type: stdout
#stdout: tibble.txt          #this will appear in the working directory
