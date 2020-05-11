library(tidyverse)
library(data.table)

#   command line
args <- commandArgs(trailingOnly = TRUE)

csv <- args[1]
bait <- args[2]
mass_tolerance <- as.numeric(args[3])
#output_path <- args[4]

#path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/ANN-Solo runs/26.11.2019/"
#name <- "b10002_TRAF1"
#mass_tolerance <- 50

#input_path <- paste(path, "csv/", name, ".csv", sep = "")
#output_path <- paste(path, "tol_py/", name, "_tol_py.csv", sep = "")
output_path <- paste(bait, "_tol_py.csv", sep = "")

#   mod file (uitkomst python) bekijken
#mod_path <- paste(path, "mod/","b10000_ZNF230_mod", ".csv", sep = "")
#mod = read.csv(mod_path, sep = '\t')
#mod %>% as_tibble %>% select(mod_mass, sequence, mass_diff, mod) %>% print(n=50)

#   mztab file
mztab = fread(file=csv, sep = '\t')

#   calculate the mass tolerance in Da for each psm
mztab_mtol <- mztab %>% as_tibble %>%
    select(sequence, PSM_ID, exp_mass_to_charge, calc_mass_to_charge, charge) %>%
    mutate(mass_diff = (exp_mass_to_charge - calc_mass_to_charge)*charge) %>%
    mutate(calc_mass = calc_mass_to_charge * charge) %>%
    mutate(mass_tol = mass_tolerance * calc_mass * 10**-6) %>%
    mutate(mass_tol_pos = mass_diff + mass_tol) %>%
    mutate(mass_tol_neg = mass_diff - mass_tol) %>%
    filter(!(mass_tol_neg < 0 & mass_tol_pos > 0)) %>%      # filter out unmodified psms
    mutate(bait = bait) %>%
    arrange(mass_tol_pos)

#   save to file for python
fwrite(mztab_mtol, output_path, na = "NA", append = FALSE, col.names = TRUE)


#######################################################################################################################################

#   unimod
#unimod_original <- read.csv("/Users/adams/Documents/master/master 2/Stage en thesis/Data/unimod/unimod.csv", 
 #   sep = ';', stringsAsFactors = FALSE)

#unimod <- unimod_original %>% as_tibble %>%
 #   select(Monoisotopic.mass, mod_des = Description, mod_name = Interim.name) %>%
  #  mutate(mod_mass = as.numeric(Monoisotopic.mass)) %>%
   # select(mod_des, mod_name, mod_mass)

#fwrite(unimod, "/Users/adams/Documents/master/master 2/Stage en thesis/Data/unimod_py.csv",
 #   na = "NA", append = FALSE, col.names = TRUE)

#######################################################################################################################################

# see which mod_masses fit in the mass tolerance boundries
#potential_mod <- mztab_mtol %>% 
 #   mutate(
  #      type = case_when(
   #     unimod$mod_mass > mass_tol_neg & unimod$mod_mass < mass_tol_pos ~ "mod",
    #    TRUE                      ~ "no mod"
     #   ))
        
       # if_else(unimod$mod_mass > mass_tol_neg & < mass_tol_pos, 
        #mod = "TRUE"))
        #select(sequence, mass_diff, exp_mass, calc_mass, unimod$mod_name, unimod$mod_dis, unimod$mod_mass),
        #select(sequence, mass_diff, exp_mass, calc_mass, mod_name="no match", mod_dis=""), missing = NULL)
