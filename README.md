Open modification searching of SARS-CoV-2–human protein interaction data reveals novel viral modification sites
========
The code in this repository has been created to post-process open modification searching results generated with the spectral library search engine ANN-SoLo <https://github.com/bittremieux/ANN-SoLo>. We reanalyzed public affinity purification mass spectrometry data using open modification searching to investigate the presence of post-translational modifications (PTMs) in the context of the SARS-CoV-2 virus–host PPI network. The preprint is available on BioRxiv <https://www.biorxiv.org/content/10.1101/2022.03.10.483652v1>.

Description of the code
-------
**CWL workflows** were created to perform protein inference (ProteinInferenceWorkflow.cwl) and to map PTMs on observed mass differences (ModificationMappingWorkflow.cwl).
They can be run by using cwltool <https://github.com/common-workflow-language/cwltool>.
```
conda install -c conda-forge nodejs
conda install -c conda-forge cwltool
```
For the protein inference workflow **pandas** and **pyopenms** need to be installed
```
conda install -c anaconda pandas
pip install pyopenms
```
You will also need to install OpenMS <https://abibuilder.informatik.uni-tuebingen.de/archive/openms/OpenMSInstaller/release/> and PIA <https://github.com/mpc-bioinformatics/pia/releases/tag/1.4.4>.

Here is an example of the arguments needed to run the protein inference workflow:
```
cwltool --outdir /proteins ProteinInferenceWorkflow.cwl --mztab_to_idxml_py Code/mztab_to_idxml.py --mztab test.mztab --filename test --fasta test.fasta --idxml_output test.idxml --cp /pia-1.3.11/pia-1.3.11.jar --compiler de.mpc.pia.intermediate.compiler.PIACompiler --idxml test.idxml --xml_output test.xml --jar /pia-1.3.11/pia-1.3.11.jar --paramFile /param/pia_pipe_SARS_v3.xml --proteinExport test.mztab --mztab_to_csv_py Code/mztab_to_protein_csv.py --pia_proteins_R Code/pia_proteins.R
```

Additionally we created some code to perform some analyses such as a GO enrichment analysis and PPI filtering. These can be found under /Analysis-code.

Contact
-------
For more information you can send an email to <charlotte.adams@uantwerpen.be>.
