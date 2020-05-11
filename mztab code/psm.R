library(tidyverse)
library(data.table)
args <- commandArgs(trailingOnly = TRUE)

csv <- args[1]
bait <- args[2]

# input_path <- paste(path, "csv/",  bait, ".csv", sep = "")
# output_path <- paste(path, "psm_csv/", bait, ".csv", sep = "")
output_path <- paste(bait, ".csv", sep = "")

csv_in = fread(csv, sep = '\t')

csv_bait <- csv_in %>% as_tibble %>%
    mutate(bait = bait)

#   save file
fwrite(csv_bait, output_path, append = FALSE, col.names = TRUE)