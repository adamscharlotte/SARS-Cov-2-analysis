# install.packages('readxl')

library(tidyverse)
library(data.table)
library(readxl)
library("cleaver")

path_table_1 <- "/Users/adams/Projects/SARS-CoV-2/Results/Supplementary/Supplementary-Table-1-v1.xlsx" # nolint

tbl_table_1 <- read_excel(path_table_1, sheet = 2)

tbl_table_cleaved <- tbl_table_1 %>%
    mutate(cleaved_sequence = cleave(sequence, enzym = "trypsin"))

 <- tbl_table_cleaved %>%
    select(PSM_ID, cleaved_sequence, everything()) %>%
    unnest(cleaved_sequence) %>%
    distinct() %>%
    add_count(PSM_ID) %>%
    mutate(missed_cleavages = n - 1) %>%
    select(-c(cleaved_sequence, n)) %>%
    distinct()

tbl_missed_cleavages %>% count(missed_cleavages)

path_csv <- "/Users/adams/Projects/SARS-CoV-2/Results/Supplementary/Supplementary-Table-1.csv" # nolint

fwrite(tbl_missed_cleavages, path_csv)