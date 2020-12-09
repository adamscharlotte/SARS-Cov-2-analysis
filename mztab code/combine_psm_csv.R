library(tidyverse)
library(data.table)

#   Load PSM data
path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/peptideIndexer/psm_cwl"
setwd(path)
files <- dir(pattern = "*.csv")
psm_csv <- files %>%
    map(read_csv) %>%       # read in all the files individually, using
                            # the function read_csv() from the readr package
    reduce(rbind)           # reduce with rbind into one dataframe

path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/csv/psm"
setwd(path)
files <- dir(pattern = "*.csv")
csv <- files %>%
    map(read_tsv) %>%       # read in all the files individually, using
                            # the function read_tsv() from the readr package
    reduce(rbind)           # reduce with rbind into one dataframe

#   Load modification data
path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/modifications"
setwd(path)
files <- dir(pattern = "*.csv")
mod_data <- files %>%
    map(read_csv) %>%       # read in all the files individually, using
                            # the function read_csv() from the readr package
    reduce(rbind) %>%       # reduce with rbind into one dataframe
    select(bait, mod, mod_mass, mass_diff, everything())

csv_map <- csv %>% select(PSM_ID, spectra_ref, exp_mass_to_charge, retention_time)
psm_map <- psm_csv %>% select(bait, sequence, accession, retention_time, exp_mass_to_charge)
psm_id <- merge(psm_map, csv_map) %>% as_tibble
psm_mod <- merge(psm_id, mod_data, by="PSM_ID") %>% as_tibble

#   Add annotation
annotation_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/names/Gordon1_annotation.txt"
annotation <- fread(annotation_path, sep= "\t")

psm_mod_annotation <- merge(psm_mod, annotation, by.x = "bait.x", by.y = "fileID") %>% 
    as_tibble %>% filter(Comments=="")        #Filter out samples that were not included in the hcip analysis.
psm_annotations <- merge(psm_id, annotation, by.x = "bait", by.y = "fileID") %>% 
    as_tibble %>% filter(Comments=="") %>%
    filter(!accession=="null")

#   Create an overview of the modifications present in the PSMs that could be matched to proteins

################################################################################################################################################
#   Only select modifications present in reported HCIPs
#   Load HCIP result file
HCIP_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/SAINT<=0.05_MiST>=0.7_AvgSpec>=2.csv"
hcip <- fread(HCIP_path) %>% as_tibble

#   Create column with BPAccession
BP_hcip <- hcip %>% 
    select(idBait, Prey, everything()) %>%
    unite(BPAccession, idBait:Prey ,sep ="_", remove=FALSE)
BP_psm_mod <- psm_mod_annotation %>% 
    separate(accession, into=c("pre","name","post","twee","prot"), remove=FALSE) %>% 
    filter(!name=="SARS") %>%
    separate(accession, into = c("pre", "Prey", "gene")) %>%
    select(Condition, Prey, everything()) %>%
    unite(BPAccession, Condition:Prey ,sep ="_", remove=FALSE) 
BP_psm <-  psm_annotations %>%
    separate(accession, into=c("pre","name","post","twee","prot"), remove=FALSE) %>% 
    filter(!name=="SARS") %>%
    separate(accession, into = c("pre", "Prey", "gene")) %>%
    select(Condition, Prey, everything()) %>%
    unite(BPAccession, Condition:Prey ,sep ="_", remove=FALSE)

psm_mod_annotation %>% 
    separate(accession, into=c("pre","name","post","twee","prot"), remove=FALSE) %>% 
    filter(name=="SARS") %>% select(mod, mass_diff, prot) %>% unique

#   Link the modification information to the hcip information
hcip_mod <- merge(BP_psm_mod, BP_hcip, by="BPAccession") %>% as_tibble
#   How many PSMs belong to a HCIP
hcip_psm <- merge(BP_psm, BP_hcip, by="BPAccession") %>% as_tibble

hcip_mod %>% filter(BP == "nsp11_COG2") %>% select(Spec, bait.x, sequence.x) %>% unique
hcip_psm %>% filter(BP == "nsp11_COG2") %>% select(Spec, bait)
hcip_mod %>% filter(BP == "nsp11_COG2") %>% select(mod, mass_diff, sequence.x)# %>% 
    #mutate(mass_diff = round(mass_diff)) %>% unique

hcip_mod %>% filter(BP == "orf9b_TOMM70") %>% select(Spec, bait.x)
hcip_psm %>% filter(BP == "orf9b_TOMM70") %>% select(Spec, bait)
hcip_mod %>% filter(BP == "orf9b_TOMM70") %>% select(mod, mass_diff, sequence.x) %>% arrange(mass_diff) %>% print(n=30)
hcip_mod %>% filter(BP == "orf9b_TOMM70") %>% count(mod) %>% arrange(desc(n)) %>% print(n=50)

hcip_psm %>% filter(BP == "orf3b_STOML2") %>% count(bait)
hcip_psm %>% filter(BP == "orf3b_STOML2", bait=="qx017172") %>% 
    select(Spec, retention_time, bait, Prey.x, PreyGene) %>% unique() %>% print(n=30)

#   There are 25 BPAccessions that do not contain any PSMs with a modifications
BP_nomod <- BP_hcip %>% filter(!BPAccession %in% hcip_mod$BPAccession)

#   Look at modifications occurring in neighbourhood of certain SARS-CoV-2 protein
hcip_mod %>% filter(Condition == "nsp8") %>% count(mod) %>% arrange(desc(n)) %>% print(n=120) 
hcip_mod %>% filter(Condition == "nsp13") %>% 
    filter(mod == "Sulfation / Phospho") %>%
    select(mod,Prey.x,PreyGene,sequence.x, PSM_ID) %>% arrange(PreyGene) %>% unique
hcip_mod %>% filter(Condition == "N") %>%
    filter(mod == "Sulfation / Phospho") %>%
    select(Prey.x,PreyGene,sequence.x, PSM_ID) %>% arrange(PreyGene) %>% unique

################################################################################################################################################

#   Look at modifications in the SARS-CoV-2 proteins
SARS_psms <- psm_annotations %>% 
    separate(accession, into=c("pre","name","post","twee","prot")) %>% 
    filter(name=="SARS") %>%
    mutate(prot=str_replace(prot,"Protein14", "orf9c(protein14)")) %>%
    filter(str_detect(Condition, prot))     #So nsp5 is matched to nsp5 and nsp5_C145A
SARS_mod <- psm_mod_annotation %>% 
    separate(accession, into=c("pre","name","post","twee","prot")) %>% 
    filter(name=="SARS") %>%
    mutate(prot=str_replace(prot,"Protein14", "orf9c(protein14)")) %>%
    filter(str_detect(Condition, prot))

#   Look at modifications per SARS-CoV-2 protein
SARS_mod %>% count(prot) %>% arrange(desc(n)) %>% print(n=30)
SARS_mod %>% filter(prot == "Protein14") %>% select(sequence.x, mod)
SARS_mod %>% filter(mod == "Sulfation / Phospho") %>% select(Condition, sequence.x, PSM_ID) %>% unique

SARS_mod %>% mutate(round_mass=round(mass_diff, 2)) %>%
    add_count(round_mass) %>% 
    select(round_mass, n, mod) %>%
    unique() %>%
    arrange(desc(n)) %>%
    # print(n=50)
    filter(n>400) %>%
    pull(round_mass)

SARS_mod %>% mutate(round_mass=round(mass_diff, 2)) %>%
    select(sequence.x, round_mass, everything()) %>%
    unite(seqmass, sequence.x:round_mass, remove=FALSE) %>%
    count(seqmass) %>% arrange(desc(n))

SARS_mod %>% count(mod) %>% 
    arrange(desc(n)) %>%
    filter(n>500) %>% pull(mod)

SARS_mod %>% select(sequence.x, mod, everything()) %>%
    unite(mod_seq, sequence.x:mod, remove=FALSE) %>%
    filter(!mod=="No direct match found in Unimod") %>%
    filter(!str_detect(mod, "Deamidation")) %>%
    count(mod_seq) %>% arrange(desc(n)) %>% 
    # filter(n>=3) %>%
    # pull(mod_seq)  
    print(n=50) 

SARS_mod %>% filter(sequence.x == "QMSCAAGTTQTACTDDNALAYYNTTK") %>% mutate(round_mass=round(mass_diff, 2)) %>%
    count(round_mass) %>% arrange(desc(n))

SARS_mod %>% filter(prot == "E") %>% #count(sequence.x) %>% arrange(desc(n)) %>% print(n=30)
    filter(sequence.x == "FNGIGVTQNVLYENQK") %>% count(mod) %>% arrange(desc(n))

SARS_mod %>% filter(prot == "nsp9") %>% #count(sequence.x) %>% arrange(desc(n)) %>% print(n=30)
    filter(sequence.x == "QMSCAAGTTQTACTDDNALAYYNTTK") %>% count(mod) %>% arrange(desc(n))

SARS_mod %>% filter(prot == "E") %>% select(sequence.y, mass_tol, mass_diff, mod) %>% pull(mod)

SARS_mod %>% as_tibble %>%
    # filter(str_detect(mod, "Gluratylation")) %>%
    select(sequence.y, mass_tol, mod, prot) %>% #count(prot)
    filter(prot == "orf3a") %>% #pull(sequence.y)
    pull(mod)
    unique %>% print(n=30)
