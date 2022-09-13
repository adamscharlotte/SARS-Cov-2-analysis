library(tidyverse)
library(data.table)


path_annsolo <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/AP-MS/ANN-SoLo/annsolo_bp.txt" # nolint
path_gordon <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/AP-MS/Gordon/gordon_bp.txt" # nolint
path_li <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/AP-MS/Li/Li_bp.txt" # nolint
path_chen <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/AP-MS/Chen/Chen_bp.txt" # nolint
path_stukalov <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/AP-MS/Stukalov/Stukalov_bp.txt" # nolint
path_germain <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/BioID/Germain/Germain_bp.txt" # nolint
path_laurent <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/BioID/Laurent/Laurent_bp.txt" # nolint
path_samavarchi <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/BioID/Samavarchi/Samavarchi_bp.txt" # nolint
path_liu <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/Liu/Liu_bp.txt" # nolint
path_bittremieux <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/AP-MS/Bittremieux/Bittremieux_bp.txt" # nolint

path_overlap <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/overlap.txt" # nolint
path_other <- "/Users/adams/Projects/SARS-CoV-2/Workspace/Other studies/other studies.txt" # nolint

# ------------------------------ Load PPI results -----------------------------

list_annsolo <- fread(path_annsolo, header = FALSE) %>%
    as_tibble() %>%
    distinct()
list_gordon <- fread(path_gordon, header = FALSE) %>%
    as_tibble() %>%
    distinct()

list_overlap <- list_annsolo %>%
    filter(V1 %in% list_gordon$V1) %>%
    separate(V1, into = c("bait", "gene"), sep = "_") %>%
    select(gene)
fwrite(list_overlap, path_overlap)

list_li <- fread(path_li, header = FALSE) %>%
    as_tibble() %>%
    filter(V1 %in% list_annsolo$V1) %>%
    distinct()
list_chen <- fread(path_gordon, header = FALSE) %>%
    as_tibble() %>%
    filter(V1 %in% list_annsolo$V1) %>%
    distinct()
list_stukalov <- fread(path_stukalov, header = FALSE) %>%
    as_tibble() %>%
    filter(V1 %in% list_annsolo$V1) %>%
    distinct()
list_germain <- fread(path_germain, header = FALSE) %>%
    as_tibble() %>%
    filter(V1 %in% list_annsolo$V1) %>%
    distinct()
list_laurent <- fread(path_laurent, header = FALSE) %>%
    as_tibble() %>%
    filter(V1 %in% list_annsolo$V1) %>%
    distinct()
list_samavarchi <- fread(path_samavarchi, header = FALSE) %>%
    as_tibble() %>%
    filter(V1 %in% list_annsolo$V1) %>%
    distinct()
list_liu <- fread(path_liu, header = FALSE) %>%
    as_tibble() %>%
    filter(V1 %in% list_annsolo$V1) %>%
    distinct()
list_bittremieux <- fread(path_bittremieux, header = FALSE) %>%
    as_tibble() %>%
    filter(V1 %in% list_annsolo$V1) %>%
    distinct()

list_other <- rbind(list_li, list_chen, list_stukalov, list_germain,
    list_laurent, list_samavarchi, list_liu, list_bittremieux) %>%
    distinct() %>%
    separate(V1, into = c("bait", "gene"), sep = "_") %>%
    select(gene)
fwrite(list_other, path_other)