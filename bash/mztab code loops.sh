#ROOT_DIR=/Users/adams/Documents/master/master\ 2/Stage\ en\ thesis/Data/ANN-Solo\ Runs

#   Add gene names
name_file=/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/mztab_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time Rscript mztab\ code/pia_genes.R "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/" "${array[0]}"
done

#   Create csv files
name_file=/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/mztab_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time python3 mztab\ code/mztab_to_psm_csv.py /Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/mztab/${array[0]}.mztab "${array[0]}"
done

name_file=/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/mztab_names.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time python3 mztab\ code/mztab_to_csv.py /Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/mztab/${array[0]}.mztab "${array[0]}"
done

#   Plot SSM

name_file=/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/SNO.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/anaconda3/envs/ann_solo/bin/python3 mztab\ code/plot_ssm.py  "${array[0]}"  "${array[1]}"  "${array[2]}"  "${array[3]}"
done

#   Plot SSM

name_file=/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/UBI.txt
lines=`tail -n+1 $name_file`
for line in $lines
do
    IFS=';' read -r -a array <<< "$line"
    time /Users/adams/anaconda3/envs/ann_solo/bin/python3 mztab\ code/plot_ssm.py  "${array[0]}"  "${array[1]}"  "${array[2]}"  "${array[3]}"
done
