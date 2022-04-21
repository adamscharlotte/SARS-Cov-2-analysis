Open modification searching of SARS-CoV-2–human protein interaction data reveals novel viral modification sites
========
The code in this repository has been created to post-process open modification searching results generated with the spectral library search engine ANN-SoLo <https://github.com/bittremieux/ANN-SoLo>. We reanalyzed public affinity purification mass spectrometry data using open modification searching to investigate the presence of post-translational modifications (PTMs) in the context of the SARS-CoV-2 virus–host PPI network. The preprint is available on BioRxiv <https://www.biorxiv.org/content/10.1101/2022.03.10.483652v1>.

Description of the code
-------
**CWL workflows** were created to perform protein inference (ProteinInferenceWorkflow.cwl) and to map PTMs on observed mass differences (ModificationMappingWorkflow.cwl).
They can be run by using cwltool <https://github.com/common-workflow-language/cwltool>. For example:
<cwltool Protein_Inference_workflow.cwl --mztab_to_idxml_py /mztab\ code/mztab_to_idxml.py --mztab --fasta uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016_sars_cov_2.fasta --idxml_output --cp pia-1.3.11.jar --compiler de.mpc.pia.intermediate.compiler.PIACompiler --bait --xml_output --jar pia-1.3.11.jar --paramFile --proteinExport --mztab_to_csv_py /mztab\ code/mztab_to_protein_csv.py --pia_proteins_R /mztab\ code/pia_proteins.R>

Additionally we created some code to perform some analyses such as a GO enrichment analysis and PPI filtering. These can be found under /Analysis-code.

Contact
-------
For more information you can send an email to <charlotte.adams@uantwerpen.be>.
