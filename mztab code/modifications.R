library(tidyverse)
library(data.table)

#   command line
args <- commandArgs(trailingOnly = TRUE)

csv <- args[1]
bait <- args[2]
mass_tolerance <- as.numeric(args[3])
output_path <- paste(bait, "_tol_py.csv", sep = "")

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