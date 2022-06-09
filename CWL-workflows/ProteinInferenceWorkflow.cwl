cwlVersion: v1.0
class: Workflow
requirements:
    MultipleInputFeatureRequirement: {}

inputs:
    mztab_to_idxml_py: File
    mztab: File
    idxml: string
    filename: string
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
    cp: File
    compiler: string
    xml_output: string
    name:
        type: string
        default: "test"
    jar: File
    proteinExport: string
    paramFile: File
    fileType: 
        type: string
        default: "mztab"    
    mztab_to_csv_py: File
    pia_proteins_R: File
    # pia_genes_R: File
    # filename: string    

outputs:
    pia_proteins:
        type: File
        outputSource: pia_proteins/proteins_csv
    # pia_genes:
    #     type: File
    #     outputSource: pia_genes/genes_csv

steps:
    mztab_to_idxml:
        run: Tools/mztab_to_idxml.cwl
        in:
            mztab_to_idxml_py: mztab_to_idxml_py
            mztab: mztab
            idxml: idxml
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

    pia_xml:
        run: Tools/pia_xml.cwl
        in:
            cp: cp
            compiler: compiler
            idxml: peptide_indexer/idxml
            filename: filename
            xml_output: xml_output
            name: name
        out:
            [xml]
    
    pia_run:
        run: Tools/pia_run.cwl
        in: 
            jar: jar
            xml: pia_xml/xml
            filename: filename
            fileType: fileType
            paramFile: paramFile
            proteinExport: proteinExport
        out:
            [mztab]
            # [csv]

    mztab_to_csv:
        run: Tools/mztab_to_csv.cwl
        in:
            mztab_to_csv_py: mztab_to_csv_py
            mztab: pia_run/mztab
            filename: filename
        out:
            [csv]

    pia_proteins:
        run: Tools/pia_proteins.cwl
        in: 
            pia_proteins_R: pia_proteins_R
            csv: mztab_to_csv/csv
            filename: filename
        out:
            [proteins_csv]

    # pia_genes:
    #     run: Tools/pia_genes.cwl
    #     in: 
    #         pia_genes_R: pia_genes_R
    #         proteins_csv: pia_proteins/proteins_csv
    #         filename: filename
    #     out:
    #         [genes_csv]
