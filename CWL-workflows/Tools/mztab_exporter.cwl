cwlVersion: v1.0
class: CommandLineTool

baseCommand: [/Applications/OpenMS-2.4.0/bin/MztabExporter]

inputs:
    input:
        type: File
        inputBinding:
            position: 1
            prefix: -in

    mztab_output:
        type: string
        inputBinding: 
            position: 2
            prefix: -out

outputs:
    mztab:
        type: File
        outputBinding: 
            glob: $(inputs.mztab_output)
