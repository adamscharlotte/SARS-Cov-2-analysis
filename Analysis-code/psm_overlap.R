library(tidyverse)
library(data.table)
library(ggvenn)

spaceless <- function(x) {
    colnames(x) <- gsub(" ", "_", colnames(x))
    x
    }

path_ann <- "/Users/adams/Projects/SARS-CoV-2/Workspace/ANN-SoLo/mztab-csv"
setwd(path_ann)
files <- dir(pattern = "*.csv")
tbl_psm_ann <- files %>%
    map(read_tsv) %>%       # read in all the files individually, using
                            # the function read_csv() from the readr package
    reduce(rbind)           # reduce with rbind into one dataframe

tbl_psm_ann_solo <- tbl_psm_ann %>%
    mutate(Spectrum = PSM_ID,
    Peptide = sequence,
    Delta_Mass = (exp_mass_to_charge - calc_mass_to_charge) * charge) %>%
    select(Spectrum, Peptide, Delta_Mass)

path_fragger <- "/Users/adams/Projects/SARS-CoV-2/Workspace/MSFragger"
tbl_psm_fragger <- fread(paste(path_fragger, "/psm.tsv", sep = "")) %>%
    as_tibble() %>%
    spaceless()

# tbl_psm_ann_solo %>%
    # filter(Spectrum %in% tbl_psm_fragger$Spectrum)
# tbl_psm_fragger %>%
#     filter(Spectrum %in% tbl_psm_ann_solo$Spectrum) %>%
#     select(Spectrum, Peptide, Delta_Mass)

path_gordon <- "/Users/adams/Projects/SARS-CoV-2/Workspace/20200318-qx-final/msms.txt" # nolint
tbl_gordon <- fread(path_gordon) %>%
    as_tibble() %>%
    spaceless()

# ----------------------------- Prepare Dataframes -----------------------------

tbl_frame_ann_solo <- tbl_psm_ann_solo %>%
    separate(Spectrum, into = c("Raw_file", "Scan_number", "scn", "charge"),
    remove = FALSE) %>%
    mutate(Scan_number = as.integer(Scan_number)) %>%
    rename(Modified_sequence = Peptide,
    Scan_ID = Spectrum) %>%
    select(Raw_file, Scan_number, Modified_sequence, Scan_ID)

tbl_frame_ann_solo$Sequence <- gsub("n.44]", "",
    gsub("n.43]", "",
    gsub(".160]", "",
    gsub(".147]", "",
    gsub(".115]", "",
    gsub(".129]", "", tbl_frame_ann_solo$Modified_sequence))))))

tbl_frame_msfragger <- tbl_psm_fragger %>%
    separate(Spectrum, into = c("Raw_file", "Scan_number", "scn", "charge"),
    remove = FALSE) %>%
    mutate(Scan_number = as.integer(Scan_number)) %>%
    rename(Sequence = Peptide,
    Scan_ID = Spectrum) %>%
    select(Raw_file, Scan_number, Sequence, Scan_ID)

tbl_frame_original <- tbl_gordon %>%
    select(Raw_file, Scan_number, Sequence, Modified_sequence)

# tbl_frame_msfragger %>%
#     arrange(Scan_number) %>%
#     arrange(Raw_file) %>%
#     print(n = 50)

psm_annsolo <- tbl_frame_ann_solo %>%
    mutate(Sequence_map = gsub("[IL]", "X", Sequence)) %>%
    select(Raw_file, Scan_number, Sequence_map, Sequence, Scan_ID) %>%
    unite(PSM, Raw_file:Sequence_map, sep = "|", remove = FALSE) %>%
    distinct()

psm_fragger <- tbl_frame_msfragger %>%
    mutate(Sequence_map = gsub("[IL]", "X", Sequence)) %>%
    select(Raw_file, Scan_number, Sequence_map, Sequence, Scan_ID) %>%
    unite(PSM, Raw_file:Sequence_map, sep = "|", remove = FALSE) %>%
    distinct()

psm_original <- tbl_frame_original %>%
    mutate(Sequence_map = gsub("[IL]", "X", Sequence)) %>%
    select(Raw_file, Scan_number, Sequence_map, Sequence) %>%
    unite(PSM, Raw_file:Sequence_map, sep = "|", remove = FALSE) %>%
    distinct()

# -------------------------------- Venn Diagram --------------------------------

list_venn <- list(
    ann_solo = psm_annsolo$PSM,
    msfragger = psm_fragger$PSM,
    original = psm_original$PSM
)

complete_venn <- c(psm_annsolo$PSM,
    psm_fragger$PSM,
    psm_original$PSM) %>%
    unique() %>%
    as_tibble()

data_venn <- data.frame(PSM = complete_venn$value,
    ann_solo = FALSE,
    msfragger = FALSE,
    original = FALSE
    )

data_venn$ann_solo <- data_venn$PSM %in% list_venn$ann_solo
data_venn$msfragger <- data_venn$PSM %in% list_venn$msfragger
data_venn$original <- data_venn$PSM %in% list_venn$original

ggplot(data_venn, aes(A = ann_solo, B = msfragger, C = original)) +
    geom_venn() +
    theme_void()

# ------------------- Plot PSM overlap MSFragger vs Original -------------------

fragger <- psm_fragger %>%
    filter(!PSM %in% psm_original$PSM) %>%
    mutate(Identification = "MSFragger") %>%
    mutate(Level = "PSM") %>%
    unique()

gordon <- psm_original %>%
    filter(!PSM %in% psm_fragger$PSM) %>%
    mutate(Scan_ID = "") %>%
    mutate(Identification = "Original analysis") %>%
    mutate(Level = "PSM") %>%
    unique()

both <- psm_fragger %>%
    filter(PSM %in% psm_original$PSM) %>%
    mutate(Identification = "Overlap") %>%
    mutate(Level = "PSM") %>%
    unique()

317164 / 388342

574394 / 388342

830743 / 388342

397308 / 574394

comb_psm <- rbind(fragger, gordon, both)
comb_psm$Identification <- factor(comb_psm$Identification, levels=c("MSFragger", "Overlap", "Original analysis")) # nolint

plot_psm <- ggplot(comb_psm, aes(Level, fill = Identification)) +
    geom_bar(position = "stack") +
    theme_minimal() +
    labs(title = "PSM overlap", x = "", y = "Number of PSMs") +
    scale_fill_manual(values = c("#C0D2F7", "#0E1C36", "#F33B16"), name = "") +
    scale_y_continuous(labels = scales::comma) +
    theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

path_plot <- "/Users/adams/Projects/SARS-CoV-2/Results/Figures/Overlap/PSM overlap - MSFragger vs Original.png" # nolint
ggsave(path_plot, plot_psm, width = 7.5, height = 8, units = "cm")

# ------------------- Plot PSM overlap ANN-SoLo vs MSFragger -------------------

ann <- psm_annsolo %>%
    filter(!PSM %in% psm_fragger$PSM) %>%
    mutate(Identification = "Original annsol") %>%
    mutate(Level = "PSM") %>%
    unique()

fragger <- psm_fragger %>%
    filter(!PSM %in% psm_annsolo$PSM) %>%
    mutate(Identification = "MSFragger") %>%
    mutate(Level = "PSM") %>%
    unique()

both <- psm_annsolo %>%
    filter(PSM %in% psm_fragger$PSM) %>%
    mutate(Identification = "Overlap") %>%
    mutate(Level = "PSM") %>%
    unique()

comb_psm <- rbind(ann, fragger, both)
comb_psm$Identification <- factor(comb_psm$Identification, levels=c("Original annsol", "Overlap", "MSFragger")) # nolint

plot_psm <- ggplot(comb_psm, aes(Level, fill = Identification)) +
    geom_bar(position = "stack") +
    theme_minimal() +
    labs(title = "PSM overlap", x = "", y = "Number of PSMs") +
    scale_fill_manual(values = c("#C0D2F7", "#0E1C36", "#F33B16"), name = "") +
    scale_y_continuous(labels = scales::comma) +
    theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

path_plot <- "/Users/adams/Projects/SARS-CoV-2/Results/Figures/Overlap/PSM overlap - MSFragger vs ANN-SoLo.png" # nolint
ggsave(path_plot, plot_psm, width = 7.5, height = 8, units = "cm")

psm_annsolo %>%
    filter(Scan_ID %in% psm_fragger$Scan_ID) %>%
    filter(!PSM %in% psm_fragger$PSM) %>%
    arrange(Scan_ID) %>%
    select(Scan_ID, PSM)

psm_fragger %>%
    filter(Scan_ID %in% psm_annsolo$Scan_ID) %>%
    filter(!PSM %in% psm_annsolo$PSM) %>%
    arrange(Scan_ID) %>%
    select(Scan_ID, PSM)

# ---------- Peptides -----------
ann <- psm_annsolo %>%
    filter(!Peptide_free %in% psm_fragger$Peptide) %>%
    select(Peptide_free) %>%
    rename(Peptide = Peptide_free) %>%
    mutate(Identification = "ANN-SoLo") %>%
    mutate(Level = "Peptide") %>%
    unique()

fragger <- psm_fragger %>%
    filter(!Peptide %in% psm_annsolo$Peptide_free) %>%
    select(Peptide) %>%
    mutate(Identification = "MSFragger") %>%
    mutate(Level = "Peptide") %>%
    unique()

both <- psm_annsolo %>%
    filter(Peptide_free %in% psm_fragger$Peptide) %>%
    select(Peptide) %>%
    mutate(Identification = "Overlap") %>%
    mutate(Level = "Peptide") %>%
    unique()

comb_pep <- rbind(ann, fragger, both)
comb_pep$Identification <- factor(comb_pep$Identification, levels=c("ANN-SoLo", "Overlap", "MSFragger")) # nolint

plot_pep <- ggplot(comb_pep, aes(Level, fill = Identification)) +
    geom_bar(position = "stack") +
    theme_minimal() +
    labs(title = "Peptide overlap", x = "", y = "Number of peptides") +
    scale_fill_manual(values = c("#C0D2F7", "#0E1C36", "#F33B16"), name = "") +
    scale_y_continuous(labels = scales::comma) +
    theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

path_plot <- "/Users/adams/Projects/SARS-CoV-2/Results/Figures/Overlap/Peptide overlap - MSFragger vs ANN-SoLo.png" # nolint
ggsave(path_plot, plot_pep, width = 7.5, height = 8, units = "cm")

# --------- Look at non-matching PSMs --------

tbl_spectra_ann <- tbl_psm_ann %>%
    mutate(Spectrum = PSM_ID,
    Peptide_ann = sequence,
    Delta_Mass_ann = (exp_mass_to_charge - calc_mass_to_charge) * charge) %>%
    select(Spectrum, Peptide_ann, Delta_Mass_ann,
        `search_engine_score[2]`, `search_engine_score[1]`)

tbl_spectra_ann$Peptide_ann_free <- gsub("n.44]", "",
    gsub("n.43]", "",
    gsub(".160]", "",
    gsub(".147]", "",
    gsub(".115]", "",
    gsub(".129]", "", tbl_spectra_ann$Peptide_ann))))))

tbl_spectra_fragger <- tbl_psm_fragger %>%
    rename(Peptide_frag = Peptide,
    Delta_Mass_frag = Delta_Mass) %>%
    select(Spectrum, Peptide_frag, Delta_Mass_frag, Hyperscore, Nextscore)

tbl_spectra_overlap <- merge(tbl_spectra_ann, tbl_spectra_fragger) %>% as_tibble() # nolint

tbl_spectra_overlap %>%
    filter(!Peptide_ann_free == Peptide_frag) %>%
    select(Peptide_ann_free, Peptide_frag, Hyperscore,
    Nextscore, `search_engine_score[2]`, `search_engine_score[1]`) %>%
    print(n = 150)

tbl_spectra_overlap %>%
    filter(!Peptide_ann_free == Peptide_frag) %>%
    arrange(`search_engine_score[2]`) %>%
    select(Peptide_ann_free, Peptide_frag, Delta_Mass_ann,
    Delta_Mass_frag) %>%
    print(n = 150)

tbl_spectra_overlap %>%
    filter(Spectrum == "qx017105.15983.15983.2") %>%
    select(Delta_Mass_ann, Delta_Mass_frag, Peptide_ann)

tbl_spectra_overlap %>%
    filter(Spectrum == "qx017181.22087.22087.2") %>%
    select(Delta_Mass_ann, Delta_Mass_frag, Peptide_ann, Peptide_frag)

tbl_spectra_overlap %>%
    filter(Spectrum == "qx017102.23076.23076.2") %>%
    select(Delta_Mass_ann, Delta_Mass_frag, Peptide_ann, Peptide_frag)


# ------ Difference in Delta Mass ---------

overlap_annsolo <- psm_annsolo %>%
    filter(PSM %in% both$PSM) %>%
    select(PSM, Delta_Mass, Peptide) %>%
    rename(Delta_Mass_ann = Delta_Mass) %>%
    distinct()

overlap_fragger <- psm_fragger %>%
    filter(PSM %in% both$PSM) %>%
    mutate(Charge = as.double(Charge)) %>%
    select(PSM, Delta_Mass) %>%
    rename(Delta_Mass_frag = Delta_Mass) %>%
    distinct()

overlap_delta <- merge(overlap_annsolo, overlap_fragger) %>% as_tibble()

overlap_delta %>% arrange(Delta_Mass_ann) %>% select(-PSM) %>% distinct() %>% print(n = 80) # nolint

tbl_overlap_counts <- overlap_delta %>%
    mutate(difference = Delta_Mass_ann - Delta_Mass_frag) %>%
    arrange(desc(difference)) %>%
    mutate(difference_round_dbl = round(difference, 2)) %>%
    mutate(difference_round = as.character(difference_round_dbl)) %>%
    count(difference_round) %>%
    distinct()

order <- tbl_overlap_counts %>%
    arrange(difference_round_dbl) %>%
    pull(difference_round_dbl) %>%
    unique() %>%
    as.character()

tbl_overlap_counts$difference_round <- factor(tbl_overlap_counts$difference_round, levels=order) # nolint

plot_overlap_difference <- ggplot(tbl_overlap_counts,
    aes(x = difference_round, y = n)) +
    geom_bar(stat = "identity") +
    geom_text(size = 3, aes(label = n), vjust = -0.4)

plot_overlap_difference

tbl_psm_ann_solo %>% select(Peptide_free) %>% print(n = 20)
# ------------------------------ Plot PSM overlap ------------------------------

psm_annsolo <- tbl_psm_ann_solo %>%
    unite(PSM, Spectrum:Peptide, sep = "|", remove = FALSE) %>%
    distinct()

psm_fragger <- tbl_psm_fragger %>%
    unite(PSM, Spectrum:Peptide, sep = "|", remove = FALSE) %>%
    distinct()

ann <- psm_annsolo %>%
    filter(!PSM %in% psm_fragger$PSM) %>%
    mutate(Identification = "ANN-SoLo") %>%
    mutate(Axis = "PSM") %>%
    unique()
fragger <- psm_fragger %>%
    filter(!PSM %in% psm_annsolo$PSM) %>%
    select(PSM, Spectrum, Peptide, Delta_Mass) %>%
    mutate(Identification = "MSFragger") %>%
    mutate(Axis = "PSM") %>%
    unique()
both <- psm_annsolo %>%
    filter(PSM %in% psm_fragger$PSM) %>%
    mutate(Identification = "Overlap") %>%
    mutate(Axis = "PSM") %>%
    unique()

comb_psm <- rbind(ann, fragger, both)

psm <- ggplot(comb_psm, aes(Axis, fill = Identification)) +
    geom_bar(position = "stack") +
    theme_minimal() +
    labs(title = "PSM overlap", x = "Raw file",
    y = "Number of identifications") +
    scale_fill_manual(values = c("#071E22", "#849B96", "#377563")) +
    theme(axis.text.x = element_text(size = 9, angle = 90))

path_plot <- "/Users/adams/Projects/SARS-CoV-2/Results/Figures/Overlap/PSM overlap - MSFragger vs ANN-SoLo.png" # nolint
ggsave(path_plot, psm, width = 31, height = 9, units = "cm")
