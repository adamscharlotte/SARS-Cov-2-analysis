#	Perform protein inference
name_file=filenames.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
	IFS=';' read -r -a array <<< "$line"
	/Users/adams/anaconda3/bin/cwltool --outdir proteins /CWL\ workflows/ProteinInferenceWorkflow.cwl --mztab_to_idxml_py /mztab\ code/mztab_to_idxml.py --mztab ${array[0]}.mztab --fasta uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016_sars_cov_2.fasta --idxml_output ${array[0]}.idxml --cp pia-1.3.11.jar --compiler de.mpc.pia.intermediate.compiler.PIACompiler --bait ${array[0]} --xml_output ${array[0]}.xml --jar /Users/adams/Downloads/pia-1.3.11/pia-1.3.11.jar --paramFile param/pia_pipe_SARS_v3.xml --proteinExport ${array[0]}.mztab --mztab_to_csv_py /mztab\ code/mztab_to_protein_csv.py --pia_proteins_R /mztab\ code/pia_proteins.R
done

#	Map mass differences to modifications in the unimod database
name_file=filenames.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
	IFS=';' read -r -a array <<< "$line"
	/Users/adams/anaconda3/bin/cwltool --outdir modifications /Users/adams/Desktop/SARS-CoV-2-analysis/CWL\ workflows/ModificationMappingWorkflow.cwl --mztab mztab/${array[0]}.mztab --bait ${array[0]} --mztab_to_csv_py /Users/adams/Desktop/SARS-CoV-2-analysis/mztab\ code/mztab_to_csv.py --modifications_R /Users/adams/Desktop/SARS-CoV-2-analysis/mztab\ code/modifications.R --modifications_faster_py /Users/adams/Desktop/SARS-CoV-2-analysis/mztab\ code/modifications_faster.py --unimod_path "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/unimod/unimod_py.csv" --mass_tolerance "10"
done

#   PSMs
name_file=filenames.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
	IFS=';' read -r -a array <<< "$line"
	/Users/adams/anaconda3/bin/cwltool --outdir peptideIndexer/psm_cwl /Users/adams/Desktop/SARS-CoV-2-analysis/CWL\ workflows//SARS_psm_workflow.cwl --mztab mztab/${array[0]}.mztab --bait ${array[0]} --mztab_to_csv_py /Users/adams/Desktop/SARS-CoV-2-analysis/mztab\ code/mztab_to_psm_csv.py --psm_R /Users/adams/Desktop/SARS-CoV-2-analysis/mztab\ code/psm.R --fasta fasta/uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016_sars_cov_2.fasta --idxml_output ${array[0]}.idxml --mztab_output ${array[0]}.mztab --mztab_to_idxml_py /Users/adams/Desktop/SARS-CoV-2-analysis/mztab\ code/mztab_to_idxml.py
done

