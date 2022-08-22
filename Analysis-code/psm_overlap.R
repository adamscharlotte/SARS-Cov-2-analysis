library(tidyverse)
library(data.table)

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

tbl_psm_ann_solo %>% filter(Spectrum %in% tbl_psm_fragger$Spectrum)
tbl_psm_fragger %>%
    filter(Spectrum %in% tbl_psm_ann_solo$Spectrum) %>%
    select(Spectrum, Peptide, Delta_Mass)

path_gordon <- "/Users/adams/Projects/SARS-CoV-2/Workspace/20200318-qx-final/msms.txt" # nolint
tbl_gordon <- fread(path_gordon) %>%
    as_tibble()


# ------------------------------ Plot PSM overlap ------------------------------

tbl_psm_ann_solo$Peptide_free <- gsub("n.44]", "",
    gsub("n.43]", "",
    gsub(".160]", "",
    gsub(".147]", "",
    gsub(".115]", "",
    gsub(".129]", "", tbl_psm_ann_solo$Peptide))))))

tbl_psm_ann_solo %>%
    filter(str_detect(Peptide_free, "]")) %>%
    select(Peptide_free) %>%
    distinct()

psm_annsolo <- tbl_psm_ann_solo %>%
    select(Spectrum, Peptide_free, everything()) %>%
    unite(PSM, Spectrum:Peptide_free, sep = "|", remove = FALSE) %>%
    distinct()

psm_fragger <- tbl_psm_fragger %>%
    select(Spectrum, Peptide, everything()) %>%
    unite(PSM, Spectrum:Peptide, sep = "|", remove = FALSE) %>%
    distinct()

ann <- psm_annsolo %>%
    filter(!PSM %in% psm_fragger$PSM) %>%
    mutate(Identification = "ANN-SoLo") %>%
    mutate(Level = "PSM") %>%
    unique()

fragger <- psm_fragger %>%
    filter(!PSM %in% psm_annsolo$PSM) %>%
    select(PSM, Spectrum, Peptide, Delta_Mass) %>%
    mutate(Identification = "MSFragger") %>%
    mutate(Level = "PSM") %>%
    unique()

both <- psm_annsolo %>%
    filter(PSM %in% psm_fragger$PSM) %>%
    mutate(Identification = "Overlap") %>%
    mutate(Level = "PSM") %>%
    unique()

comb_psm <- rbind(ann, fragger, both)
comb_psm$Identification <- factor(comb_psm$Identification, levels=c("ANN-SoLo", "Overlap", "MSFragger")) # nolint

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

# ------------- Spectra ----------

ann <- psm_annsolo %>%
    filter(!Spectrum %in% psm_fragger$Spectrum) %>%
    select(Spectrum, Peptide_free, Delta_Mass) %>%
    rename(Peptide = Peptide_free) %>%
    mutate(Identification = "ANN-SoLo") %>%
    mutate(Level = "Spectrum") %>%
    unique()

fragger <- psm_fragger %>%
    filter(!Spectrum %in% psm_annsolo$Spectrum) %>%
    select(Spectrum, Peptide, Delta_Mass) %>%
    mutate(Identification = "MSFragger") %>%
    mutate(Level = "Spectrum") %>%
    unique()

both <- psm_annsolo %>%
    filter(Spectrum %in% psm_fragger$Spectrum) %>%
    select(Spectrum, Peptide_free, Delta_Mass) %>%
    rename(Peptide = Peptide_free) %>%
    mutate(Identification = "Overlap") %>%
    mutate(Level = "Spectrum") %>%
    unique()

comb_spec <- rbind(ann, fragger, both)
comb_spec$Identification <- factor(comb_spec$Identification, levels=c("ANN-SoLo", "Overlap", "MSFragger")) # nolint

plot_spec <- ggplot(comb_spec, aes(Level, fill = Identification)) +
    geom_bar(position = "stack") +
    theme_minimal() +
    labs(title = "Spectrum overlap", x = "", y = "Number of spectra") +
    scale_fill_manual(values = c("#C0D2F7", "#0E1C36", "#F33B16"), name = "") +
    scale_y_continuous(labels = scales::comma) +
    theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

path_plot <- "/Users/adams/Projects/SARS-CoV-2/Results/Figures/Overlap/Spectrum overlap - MSFragger vs ANN-SoLo.png" # nolint
ggsave(path_plot, plot_spec, width = 7.5, height = 8, units = "cm")

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
