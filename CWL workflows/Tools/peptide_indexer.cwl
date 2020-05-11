cwlVersion: v1.0
class: CommandLineTool

baseCommand: [/Applications/OpenMS-2.4.0/bin/PeptideIndexer]

inputs:
    input:
        type: File
        inputBinding:
            position: 1
            prefix: -in
    
    fasta:
        type: File
        inputBinding:
            position: 2
            prefix: -fasta

    idxml_output:
        type: string
        inputBinding:
            position: 3
            prefix: -out

    missing_decoy_action:
        type: string
        inputBinding:
            position: 4
            prefix: -missing_decoy_action

    IL_equivalent:
        type: boolean
        inputBinding:
            position: 5
            prefix: -IL_equivalent

    enzyme_specificity:
        type: string
        inputBinding:
            position: 6
            prefix: -enzyme:specificity
    
    allow_unmatched:
        type: boolean
        inputBinding:
            position: 7
            prefix: -allow_unmatched

outputs:
    idxml:
        type: File
        outputBinding: 
            glob: $(inputs.idxml_output)
