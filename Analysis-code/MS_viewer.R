library(tidyverse)
library(data.table)

path_mztab <- "/Users/adams/Projects/SARS-CoV-2/Workspace/ANN-SoLo/mztab-csv"
setwd(path_mztab)
files <- dir(pattern = "*.csv")
tbl_mztab <- files %>%
    map(read_tsv) %>%       # read in all the files individually, using
                            # the function read_csv() from the readr package
    reduce(rbind)           # reduce with rbind into one dataframe

tbl_result <- tbl_mztab %>%
    separate(PSM_ID, into = c("file_name", "scan", "b", "c"),
    remove = FALSE) %>%
    select(-c(b, c)) %>%
    mutate(fraction = file_name) %>%
    mutate(Spec_ID = paste("controllerType=0 controllerNumber=1 scan=", scan,
    sep = "")) %>%
    select(fraction, Spec_ID, sequence, charge, everything())

colnames(tbl_result)

tbl_result$modified_sequence <- gsub("n.44]", "Carbamyl-",
    gsub("n.43]", "Acetyl-",
    gsub(".160]", "(Carbamidomethyl)",
    gsub(".147]", "(Oxidation)",
    gsub(".115]", "(Deamidated)",
    gsub(".129]", "(Deamidated)", tbl_result$sequence))))))

tbl_result <- tbl_result %>%
    select(fraction, Spec_ID, modified_sequence, charge, everything())

path_ms_viewer <- "/Users/adams/Projects/SARS-CoV-2/Workspace/ANN-SoLo/MS-Viewer.csv" # nolint
fwrite(tbl_result, path_ms_viewer)

# --------------------------------- MSFragger ----------------------------------

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

path_fragger <- "/Users/adams/Projects/SARS-CoV-2/Workspace/MSFragger"
tbl_psm_fragger <- fread(paste(path_fragger, "/psm.tsv", sep = "")) %>%
    as_tibble() %>%
    spaceless()

tbl_psm_fragger %>% select(Assigned_Modifications) %>% distinct()