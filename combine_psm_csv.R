library(tidyverse)
library(data.table)

path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/peptideIndexer/psm_cwl"
setwd(path)
files <- dir(pattern = "*.csv")
psm_csv <- files %>%
    map(read_csv) %>%       # read in all the files individually, using
                            # the function read_tsv() from the readr package
    reduce(rbind)           # reduce with rbind into one dataframe

path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/csv"
setwd(path)
files <- dir(pattern = "*.csv")
csv <- files %>%
    map(read_tsv) %>%       # read in all the files individually, using
                            # the function read_tsv() from the readr package
    reduce(rbind)           # reduce with rbind into one dataframe

csv_map <- csv %>% select(PSM_ID, spectra_ref, exp_mass_to_charge, retention_time)

SARS_psms <- psm_csv %>% separate(accession, into=c("pre","name","post","twee", "prot")) %>% 
    filter(name=="SARS") %>%
    select(bait, sequence, prot, retention_time, exp_mass_to_charge)

SARS_map <- merge(SARS_psms, csv_map) %>% as_tibble

#   modification data
path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2//modifications_20"
setwd(path)
files <- dir(pattern = "*.csv")
pre_data <- files %>%
    map(read_csv) %>%       # read in all the files individually, using
                            # the function read_tsv() from the readr package
    reduce(rbind) %>%       # reduce with rbind into one dataframe
    select(bait, mod, mod_mass, mass_diff, everything())

data_psms <- merge(SARS_map, pre_data, by="PSM_ID") %>% as_tibble

data_psms %>% select(sequence.x, sequence.y)
data_psms %>% count(prot) %>% arrange(desc(n)) %>% print(n=30)
data_psms %>% mutate(round_mass=round(mass_diff, 2)) %>%
    add_count(round_mass) %>% 
    select(round_mass, n, mod) %>%
    unique() %>%
    arrange(desc(n)) %>%
    # print(n=50)
    filter(n>400) %>%
    pull(round_mass)

data_psms %>% mutate(round_mass=round(mass_diff, 2)) %>%
    select(sequence.x, round_mass, everything()) %>%
    unite(seqmass, sequence.x:round_mass, remove=FALSE) %>%
    count(seqmass) %>% arrange(desc(n))

data_psms %>% count(mod) %>% 
    arrange(desc(n)) %>%
    filter(n>500) %>% pull(mod)


count(sequence.x) %>% arrange(desc(n))


data_psms %>% select(sequence.x, mod, everything()) %>%
    unite(mod_seq, sequence.x:mod, remove=FALSE) %>%
    filter(!mod=="No direct match found in Unimod") %>%
    filter(!str_detect(mod, "Deamidation")) %>%
    count(mod_seq) %>% arrange(desc(n)) %>% 
    # filter(n>=3) %>%
    # pull(mod_seq)  
    print(n=50) 

data_psms %>% filter(sequence.x == "QMSCAAGTTQTACTDDNALAYYNTTK") %>% mutate(round_mass=round(mass_diff, 2)) %>%
    count(round_mass) %>% arrange(desc(n))

data_psms %>% filter(prot == "E") %>% #count(sequence.x) %>% arrange(desc(n)) %>% print(n=30)
    filter(sequence.x == "FNGIGVTQNVLYENQK") %>% count(mod) %>% arrange(desc(n))

data_psms %>% filter(prot == "nsp9") %>% #count(sequence.x) %>% arrange(desc(n)) %>% print(n=30)
    filter(sequence.x == "QMSCAAGTTQTACTDDNALAYYNTTK") %>% count(mod) %>% arrange(desc(n))

data_psms %>% filter(prot == "E") %>% select(sequence.y, mass_tol, mass_diff, mod) %>% pull(mod)

data_psms %>% as_tibble %>%
    #filter(str_detect(mod, "Gluratylation")) %>%
    select(sequence.y, mass_tol, mass_diff, mod, prot) %>% #count(prot)
    filter(prot == "Protein14") %>% pull(sequence.y)

################################################################################################################################################

data_psms %>% as_tibble %>% add_count(mod) %>% arrange(desc(n)) %>% select(mod,n) %>% 
    # unique %>% filter(n>10) %>% pull(mod)
    # filter(!str_detect(mod, "No direct match found in Unimod")) %>%
    # filter(!str_detect(mod, "Deamidation")) %>%
    # filter(!str_detect(mod, "Pyro_glu")) %>% 
    # filter(!str_detect(mod, "Acetyl")) %>%
    # filter(!str_detect(mod, "Hydroxylation")) %>%
    # filter(!str_detect(mod, "Formyl")) %>%
    # filter(!str_detect(mod, "Amide")) %>%
    # filter(!str_detect(mod, "Pyro-glu")) %>%
    # filter(!str_detect(mod, "Formyl")) %>%
    filter(str_detect(mod, "Ammonium")) %>%
    unique() %>% filter(n > 1000) %>% print(n = 40) %>% pull(mod) 

#   Histogram
ylabel <- paste("Number of PSMs")
mass_diff_histogram <- ggplot(data_psms, aes(mass_diff)) + geom_histogram(binwidth = 1) + 
    theme_minimal() +
    xlim(-60,120) +
    xlab('Precursor mass difference (Da)') + 
    ylab(ylabel)

ggsave("/Users/adams/Documents/master/master 2/Stage en thesis/Reports/Figures/mass difference histogram.png", mass_diff_histogram)


################################################################################################################################################

#   load identification data
path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/identifications"
setwd(path)
files <- dir(pattern = "*.csv")

ANN_SoLo_data <- files %>%
    map(read_csv) %>%       # read in all the files individually, using
                            # the function read_csv() from the readr package
    reduce(rbind) %>%       # reduce with rbind into one dataframe
    mutate(peptide = ifelse(modified_peptide == FALSE, "is not modified", "is modified")) %>%
    mutate(iden="ann")

iden_data <- ANN_SoLo_data %>% separate(PSM_ID, into=c("b", "PSM_ID", "id2", "num")) %>%
    select(bait, PSM_ID, everything()) %>%
    unite(bait_ID, bait:PSM_ID, remove=FALSE)

iden_data_psms <- merge(iden_data, SARS_psms, by="bait_ID") %>% as_tibble

SARS_psms

iden_data_psms %>% select(sequence.x, sequence.y)

data_psms %>% select(sequence.x, sequence.y)

