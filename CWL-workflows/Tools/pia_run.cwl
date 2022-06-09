cwlVersion: v1.0
class: CommandLineTool

baseCommand: [java]

requirements:
 - class: InlineJavascriptRequirement

inputs:
    jar:
        type: File
        inputBinding:
            position: 1
            prefix: -jar

    xml:
        type: File
        inputBinding:
            position: 2
            prefix: -infile
    
    paramFile:
        type: File
        inputBinding:
            position: 3
            prefix: -paramFile
        
    proteinExport:
        type: string
        inputBinding:
            position: 4
            prefix: -proteinExport
    
    fileType:
        type: string
        inputBinding:
            position: 5

    filename:
        type: string
        inputBinding:
            position: 6

outputs:
    mztab:
        type: File
        outputBinding: 
            glob: $(inputs.filename).mztab
    # csv:
    #     type: File
    #     outputBinding: 
    #         glob: $(inputs.filename).csv
#    tibbles:
 #       type: stdout
#stdout: tibble.txt          #this will appear in the working directory
