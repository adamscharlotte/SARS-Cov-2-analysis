cwlVersion: v1.0
class: Workflow
requirements:
    MultipleInputFeatureRequirement: {}

inputs:
    mztab_to_csv_py: File
    modifications_R: File
    mztab: File
    bait: string
    mass_tolerance:
        type: string
        default: "20"
    modifications_faster_py: File
    unimod_path: string

outputs:
    mod_csv:
        type: File
        outputSource: modifications_faster/mod_csv
        
steps:
    mztab_to_csv:
        run: Tools/mztab_to_csv.cwl
        in:
            mztab_to_csv_py: mztab_to_csv_py
            mztab: mztab
            bait: bait
        out:
            [csv]

    modifications:
        run: Tools/modifications.cwl
        in: 
            modifications_R: modifications_R
            csv: mztab_to_csv/csv
            bait: bait
            mass_tolerance: mass_tolerance
        out:
            [tol_py_csv]

    modifications_faster:
        run: Tools/modifications_faster.cwl
        in: 
            modifications_faster_py: modifications_faster_py
            bait: bait
            unimod_path: unimod_path
            tol_py_csv: modifications/tol_py_csv
        out:
            [mod_csv]