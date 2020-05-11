library(tidyverse)

#args <- commandArgs(trailingOnly = TRUE)
#path <- args[1]

# path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/ANN-Solo runs/preliminary/modifications"
path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2//modifications"

setwd(path)
files <- dir(pattern = "*.csv")

pre_data <- files %>%
    map(read_csv) %>%       # read in all the files individually, using
                            # the function read_tsv() from the readr package
    reduce(rbind) %>%       # reduce with rbind into one dataframe
    select(bait, mod, mod_mass, mass_diff, everything())

pre_data %>% as_tibble %>% add_count(mod) %>% arrange(desc(n)) %>% select(mod, n) %>% unique %>% print(n=30)


#   look at most common modifications
pre_data %>% as_tibble %>% add_count(mod) %>% arrange(desc(n)) %>% select(mod,n) %>% 
    # unique %>% filter(n>500) %>% pull(mod)
    # filter(!str_detect(mod, "No direct match found in Unimod")) %>%
    # filter(!str_detect(mod, "Deamidation")) %>%
    # filter(!str_detect(mod, "Pyro_glu")) %>% 
    # filter(!str_detect(mod, "Acetyl")) %>%
    # filter(!str_detect(mod, "Hydroxylation")) %>%
    # filter(!str_detect(mod, "Formyl")) %>%
    # filter(!str_detect(mod, "Amide")) %>%
    # filter(!str_detect(mod, "Pyro-glu")) %>%
    # filter(!str_detect(mod, "Formyl")) %>%
    # filter(!str_detect(mod, "Ammonium")) %>%
    unique() %>% filter(n > 1000) %>% print(n = 40) %>% pull(mod) 

#   create bar plot with most prevalent modifications
mod_1 <- pre_data %>% as_tibble %>% filter(!str_detect(mod, "Deamidation_O18")) %>% filter(str_detect(mod, "Deamidation")) %>% mutate(modification = "Deamidation")
mod_2 <- pre_data %>% as_tibble %>% filter(str_detect(mod, "Pyro_glu")) %>% mutate(modification = "Pyro-glu from E")
mod_3 <- pre_data %>% as_tibble %>% filter(str_detect(mod, "Pyro-glu")) %>% mutate(modification = "Pyro-glu from Q")
mod_4 <- pre_data %>% as_tibble %>% filter(!str_detect(mod, "Formylasparagine")) %>% filter(str_detect(mod, "Formyl")) %>% mutate(modification = "Formylation")
mod_5 <- pre_data %>% as_tibble %>% filter(str_detect(mod, "Hydroxylation")) %>% mutate(modification = "Hydroxylation")

# mod_5 %>% select(mod)%>% unique

mod_comb <- rbind(mod_1, mod_2, mod_3, mod_4,mod_5)
mod_comb$modification <- factor(mod_comb$modification,levels = c("Deamidation", "Formylation", "Pyro-glu from E", "Hydroxylation","Pyro-glu from Q"))

s <- ggplot(mod_comb, aes(modification)) +
    geom_bar(position = "stack", fill="#071E22") +
    theme_minimal() +
    labs(title ="Top 5 most frequent PTMs", x = "modification", y = "number of identifications") + 
    #scale_fill_manual(values=c("#849B96", "#849B96",  "#849B96", "#849B96", "#849B96")) +
    theme(axis.text.x = element_text(size=10, angle=90))

ggsave("/Users/adams/Documents/master/master 2/Stage en thesis/Reports/Figures/Top 5 most frequent PTMs.png", s, width = 8, height = 15, units = "cm")


################################################################################################################################################

#   specifically look at deamidation
num <- pre_data %>% as_tibble %>% filter(str_detect(mod, "Lys->Glu / Xle->Asn / Deamidation / Asn->Asp")) %>% 
    #select(mod) %>% count(mod)
    pull(mass_diff)
median(num)

#   specifically look at methylation
pre_data %>% as_tibble %>% filter(str_detect(mod, "Methylation")) %>% 
    select(mod) %>% count(mod)

#   specifically look at pyro_glu
num <- pre_data %>% as_tibble %>% filter(str_detect(mod, "Pyro_glu")) %>% 
    #select(mod, mod_mass) %>% add_count(mod) %>% unique() %>% pull(n)
    pull(mass_diff)
median(num)

#   specifically look at Hydroxylation
num <- pre_data %>% as_tibble %>% filter(str_detect(mod, "Hydroxylation / Ala->Ser / Phe->Tyr ")) %>% 
    #select(mod, mod_mass) %>% count(mod) %>% unique() %>% pull(n) 
    pull(mass_diff)
median(num)

#   specifically look at Amide
num <- pre_data %>% as_tibble %>% filter(str_detect(mod, "Lysaminoadipicsealde / Cystine / Amide")) %>% 
    #select(mod, mod_mass) %>% add_count(mod) %>% unique() %>% pull(n) 
    pull(mass_diff)
median(num)

#   specifically look at Amide
num <- pre_data %>% as_tibble %>% filter(str_detect(mod, "semialdehyde")) %>% 
    select(mod, mod_mass) %>% add_count(mod) %>% unique() %>% pull(n)
sum(num)

 
#   extract [...] within the sequences
pre_data %>% as_tibble %>%
    separate(sequence, into = c("pre", "mid", "post")) %>%
    select(mid) %>% count(mid) %>% arrange(desc(n))

#   specifically look at Acetylation
#   select 43 in sequence
pre_data %>% as_tibble %>% filter(str_detect(sequence, "43")) %>% 
    select(mod) %>% count(mod) %>% arrange(desc(n))

num <- pre_data %>% as_tibble %>% filter(str_detect(mod, "Acetyl / Ser->Glu")) %>% 
    #filter(!str_detect(sequence, "43")) %>%
    #select(mod) %>% count(mod) %>% arrange(desc(n)) %>% print(n=50)   #pull(mod)
    #select(mod, mod_mass) %>% add_count(mod) %>% arrange(desc(n)) %>% unique()
    pull(mass_diff)
median(num)

#   specifically look at Formylation
pre_data %>% as_tibble %>% filter(str_detect(mod, "Formyl")) %>% 
    select(mod, mod_mass) %>% add_count(mod) %>% arrange(desc(n)) %>% unique() %>% pull(mod)

#   specifically look at Dehydration
pre_data %>% as_tibble %>% filter(str_detect(mod, "Hydrox")) %>% 
    #filter(!str_detect(sequence, "43")) %>%
    #select(mod) %>% count(mod) %>% arrange(desc(n)) %>% print(n=50)   #pull(mod)
    select(mod, mod_mass) %>% add_count(mod) %>% arrange(desc(n)) %>% unique() %>% pull(mod)

#   specifically look at SNO
pre_data %>% as_tibble %>% filter(str_detect(mod, "SNO")) %>% 
    select(sequence) %>% filter(str_detect(sequence, "C")) %>% print(n=100)

pre_data %>% as_tibble %>% filter(str_detect(mod, "didehydro")) %>% 
    select(sequence) %>% filter(str_detect(sequence, "V")) %>% print(n=100)

#   add 'no direct match found'
data <- pre_data %>% as_tibble %>% 
    replace_na(list(mod = "No direct match found in Unimod", mod_mass = " "))
    
data %>% filter(mass_diff > 101.0) %>% arrange(mass_diff)
data %>% filter(mod %in% "EQIGG")


#   show frequency
data %>% select(mod, mod_mass) %>% add_count(mod) %>% arrange(desc(n)) %>% unique() %>% 
    filter(n>2200) %>% pull(mod_mass)
 
    #select(mod_mass) %>% print(n = 30)

#   select mass_diff of which "No direct match found in Unimod"
no_match_data <- data %>% as_tibble %>% filter(mod == "No direct match found in Unimod") %>% mutate(rounded = round(mass_diff, 4))
mass <- no_match_data %>% add_count(rounded, sort = TRUE) %>% filter(n > 30) %>% pull(rounded) %>% unique()
seq <- no_match_data %>% add_count(rounded, sort = TRUE) %>% filter(n > 30) %>% pull(sequence) %>% unique()

no_match_data %>% filter(rounded %in% mass) %>% 
    add_count(rounded, sort = TRUE) %>% 
    mutate(mass_tol_rounded = round(mass_tol, 4)) %>%
    mutate(calc = calc_mass - charge) %>%
    mutate(exp = exp_mass - charge) %>%
    #select(rounded, calc_mass, exp_mass, sequence, n, mass_tol_rounded, calc) %>% unique() %>% 
    select(rounded, sequence) %>% unique() %>% 
    pull(rounded)

data %>% as_tibble %>% filter(sequence %in% seq) %>% #filter(mod != "No direct match found in Unimod") %>%
    mutate(rounded = round(mass_diff, 3)) %>% add_count(rounded, sort = TRUE) %>% filter(n > 10) %>%
    select(sequence, rounded, n, mod, mod_mass, mass_tol) %>% unique() %>% print(n=50) %>% arrange(sequence) %>%
    filter(sequence  == "SNMFPNHYR") %>%
    pull(rounded)

no_match_data <- data %>% as_tibble %>% filter(mod == "No direct match found in Unimod")
#   Histogram no match data
ylabel <- paste("Number of PSMs")

no_match_histogram <- ggplot(no_match_data, aes(mass_diff)) + geom_histogram(binwidth = 1) + 
    theme_minimal() +
    xlab('Precursor mass difference (Da)') + 
    ylab(ylabel)

ggsave("/Users/adams/Documents/master/master 2/Stage en thesis/Reports/Figures/no match histogram.png", mass_diff_histogram)


#   Histogram
ylabel <- paste("Number of PSMs")

mass_diff_histogram <- ggplot(data, aes(mass_diff)) + geom_histogram(binwidth = 1) + 
    theme_minimal() +
    xlim(-60,120) +
    xlab('Precursor mass difference (Da)') + 
    ylab(ylabel)

ggsave("/Users/adams/Documents/master/master 2/Stage en thesis/Reports/Figures/mass difference histogram.png", mass_diff_histogram)



s <- ggplot(data, aes(name, fill = peptide)) +
    geom_bar(position = "stack") +
    theme_minimal() +
    labs(title ="modified vs unmodified overview", x = "file names", y = "number of identifications") + 
    scale_fill_manual(values=c("red", "pink")) +
    theme(axis.text.x = element_text(size=7, angle=90))

#   save plot
ggsave("/Users/adams/Documents/master/master 2/Stage en thesis/Reports/Figures/modified vs unmodified overview.png", s)


#   add ppm of the difference between mod mass and the obseverved mass difference
data %>% as_tibble %>%
    mutate(mod_ppm = ((as.numeric(mass_diff) - as.numeric(mod_mass))/calc_mass) *10^6) %>%
    select(mod_ppm, everything()) %>%
    print(n=50)