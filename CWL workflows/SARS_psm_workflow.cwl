cwlVersion: v1.0
class: Workflow
requirements:
    MultipleInputFeatureRequirement: {}

inputs:
    mztab_to_idxml_py: File
    mztab: File
    bait: string
    fasta: File
    idxml_output: string
    missing_decoy_action:
        type: string
        default: "silent"
    enzyme_specificity:
        type: string
        default: "full"
    IL_equivalent:
        type: boolean
        default: TRUE
    allow_unmatched:
        type: boolean
        default: TRUE 
    mztab_output: string
    mztab_to_csv_py: File
    psm_R: File

outputs:
    psm:
        type: File
        outputSource: psm/csv    

steps:
    mztab_to_idxml:
        run: Tools/mztab_to_idxml.cwl
        in:
            mztab_to_idxml_py: mztab_to_idxml_py
            mztab: mztab
            bait: bait
        out:
            [idxml]

    peptide_indexer:
        run: Tools/peptide_indexer.cwl
        in:
            input: mztab_to_idxml/idxml
            fasta: fasta
            idxml_output: idxml_output
            missing_decoy_action: missing_decoy_action
            IL_equivalent: IL_equivalent
            enzyme_specificity: enzyme_specificity
            allow_unmatched: allow_unmatched
        out:
            [idxml]

    mztab_exporter:
        run: Tools/mztab_exporter.cwl
        in:
            input: peptide_indexer/idxml
            mztab_output: mztab_output
        out:
            [mztab]

    mztab_to_csv:
        run: Tools/mztab_to_csv.cwl
        in: 
            mztab_to_csv_py: mztab_to_csv_py
            mztab: mztab_exporter/mztab
            bait: bait
        out:
            [csv]

    psm:
        run: Tools/psm.cwl
        in:
            psm_R: psm_R
            csv: mztab_to_csv/csv
            bait: bait
        out:
            [csv]
