library(tidyverse)
library(data.table)

#   Gene names
fasta_headers <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/fasta_headers.csv"
unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble
gene_names <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% drop_na(gene) %>%
    separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
    select(V2, gene_name)

#   output MiST
# Before you can load the MistOutput the first couple of lines must be removed.
MIST <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/MiST/MistOutput.txt", sep= "\t") %>% as_tibble

MIST_result <- MIST %>% filter(MiST >= 0.7) %>%
    mutate(idBait=Bait) %>%
    select(idBait, Prey, everything()) %>%
    unite(BP, idBait:Prey ,sep ="_", remove=FALSE)

mist_gene <- merge(MIST, gene_names, by.x = "Prey", by.y = "V2") %>% as_tibble

MIST_gene <- mist_gene %>% filter(MiST >= 0.7) %>%
    mutate(idBait=Bait) %>%
    mutate(PreyGene = gene_name) %>%
    select(idBait, PreyGene, everything()) %>%
    unite(BP, idBait:PreyGene ,sep ="_", remove=FALSE)

#   output SAINT
SAINT <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/SAINT/list.txt", sep= "\t") %>% as_tibble

SAINT_result <- SAINT %>% filter(BFDR <= 0.05 & AvgSpec >= 2) %>%
    mutate(idBait=Bait) %>%
    select(idBait, Prey, everything()) %>%
    unite(BP, idBait:Prey ,sep ="_", remove=FALSE)

SAINT_gene <- SAINT %>% filter(BFDR <= 0.05 & AvgSpec >= 2) %>%
    mutate(idBait=Bait) %>%
    select(idBait, PreyGene, everything()) %>%
    unite(BP, idBait:PreyGene ,sep ="_", remove=FALSE)

#   combine MIST and SAINT output
combined_result <- SAINT_result %>% filter(BP %in% MIST_result$BP)      #in both directions same amount
combined_gene <- SAINT_gene %>% filter(BP %in% MIST_gene$BP)            #in both directions same amount

output_path <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Spreadsheets/HCIP/SAINT<=0.05_MiST>=0.7_AvgSpec>=2.csv"
fwrite(combined_gene, output_path, append = FALSE, col.names = TRUE) 

################################################################################################################################################

#   load original Gordon et al. results
intact <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Gordon/intact_gordon.txt", sep= "\t")
intact_map <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Gordon/name_map_sars_cov_2.txt", sep= "\t")
intact_annot <- merge(intact, intact_map, by.x = "IDinteractorA", by.y = "accession") %>% as_tibble

#   remove space in protein14 (bait)
intact_annot$idBait <- gsub(" ", "", intact_annot$idBait)

intact_result <- intact_annot %>%
    select(idBait, IDinteractorB, everything()) %>%
    unite(BP, idBait:IDinteractorB,sep ="_", remove=FALSE) #%>% select(BP) %>% unique()

intact_G <- merge(intact_annot, gene_names, by.x = "IDinteractorB", by.y = "V2") %>% as_tibble
intact_gene <- intact_G %>%
    select(idBait, gene_name, everything()) %>%
    unite(BP, idBait:gene_name,sep ="_", remove=FALSE)

#   create plot with overlap
ann <- combined_gene %>% filter(!BP %in% intact_gene$BP) %>% select(idBait, BP) %>% mutate(Identification = "ANN-SoLo")
gordon <- intact_gene %>% filter(!BP %in% combined_gene$BP) %>% select(idBait, BP) %>% mutate(Identification = "Gordon et al.")
both <- combined_gene %>% filter(BP %in% intact_gene$BP) %>% select(idBait, BP) %>% mutate(Identification = "Overlap")

comb <- rbind(ann, gordon, both)

s <- ggplot(comb, aes(idBait, fill = Identification)) +
    geom_bar(position = "stack") +
    theme_minimal() +
    labs(title ="HCIP overlap", x = "Baits", y = "Number of identifications") + 
    scale_fill_manual(values=c("#071E22", "#849B96", "#377563")) +
    theme(axis.text.x = element_text(size=10, angle=90))

ggsave("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Results/Figures/HCIP/HCIP overlap.png", s,width = 24, height = 9, units = "cm")

# #   create plot with overlap accession
# ann <- combined_result %>% filter(!BP %in% intact_result$BP) %>% mutate(identification = "ANN-SoLo") %>% select(idBait, BP, identification)
# gordon <- intact_result %>% filter(!BP %in% combined_result$BP) %>% mutate(identification = "Gordon et al.") %>% select(idBait, BP, identification)
# both <- combined_result %>% filter(BP %in% intact_result$BP) %>% mutate(identification = "Overlap") %>% select(idBait, BP, identification)

# comb <- rbind(ann, gordon, both)

# s <- ggplot(comb, aes(idBait, fill = identification)) +
#     geom_bar(position = "stack") +
#     theme_minimal() +
#     labs(title ="HCIP overlap", x = "baits", y = "number of identifications") + 
#     scale_fill_manual(values=c("#071E22", "#849B96", "#377563")) +
#     theme(axis.text.x = element_text(size=10, angle=90))

# ggsave("/Users/adams/Documents/master/master 2/Stage en thesis/Reports/Figures/HCIP overlap.png", s,width = 22, height = 9, units = "cm")

################################################################################################################################################
combined_nsp5 <- combined_gene %>% filter(idBait == "nsp5") 
intact_nsp5 <- intact_gene %>% filter(idBait == "nsp5") %>%
    rename(PreyGene = gene_name) %>%
    select(idBait, BP, PreyGene)
all_nsp5 <- rbind(combined_nsp5, intact_nsp5) %>% unique

combined_nsp5_C145A <- combined_gene %>% filter(idBait == "nsp5_C145A")
intact_nsp5_C145A <- intact_gene %>% filter(idBait == "nsp5_C145A") %>%
    rename(PreyGene = gene_name) %>%
    select(idBait, BP, PreyGene)
all_nsp5_C145A <- rbind(combined_nsp5_C145A, intact_nsp5_C145A) %>% unique

combined_N <- combined_gene %>% filter(idBait == "N") 
intact_N <- intact_gene %>% filter(idBait == "N") %>%
    rename(PreyGene = gene_name) %>%
    select(idBait, BP, PreyGene)
all_N <- rbind(combined_N, intact_N) %>% unique

combined_nsp14 <- combined_gene %>% filter(idBait == "nsp14") 
intact_nsp14 <- intact_gene %>% filter(idBait == "nsp14") %>%
    rename(PreyGene = gene_name) %>%
    select(idBait, BP, PreyGene)
all_nsp14 <- rbind(combined_nsp14,intact_nsp14) %>% unique

output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/nsp5.csv"
fwrite(all_nsp5, output_path, append = FALSE, col.names = TRUE) 
output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/nsp5_C145A.csv"
fwrite(all_nsp5_C145A, output_path, append = FALSE, col.names = TRUE) 
output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/N.csv"
fwrite(all_N, output_path, append = FALSE, col.names = TRUE) 
output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/nsp14.csv"
fwrite(all_nsp14, output_path, append = FALSE, col.names = TRUE) 

combined_gene %>% filter(idBait == "nsp14") 
both %>% filter(idBait == "nsp14") %>% select(BP)
gordon %>% filter(idBait == "nsp14") %>% select(BP)

################################################################################################################################################

#   reproduce
nrow(both)/nrow(intact_result)
nrow(both)/nrow(intact_gene)

nrow(combined_result)/nrow(intact_result) -1
nrow(combined_gene)/nrow(intact_gene) -1

################################################################################################################################################
#   write gene names to files, to color the nodes in network

ann <- combined_gene %>% filter(!BP %in% intact_gene$BP) %>% mutate(identification = "ANN-SoLo") %>% select(idBait, BP, identification, PreyGene)
gordon <- intact_gene %>% filter(!BP %in% combined_gene$BP) %>% mutate(identification = "Gordon et al.") %>% 
    mutate(PreyGene = gene_name) %>%
    select(idBait, BP, PreyGene)
both <- combined_gene %>% filter(BP %in% intact_gene$BP) %>% mutate(identification = "Overlap") %>% select(idBait, BP, identification, PreyGene)

both_gene <- both %>% select(PreyGene) %>% unique()
ann_gene <- ann %>% select(PreyGene) %>% unique()
gordon_gene <- gordon %>% select(PreyGene) %>% unique()

output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/both_gene.csv"
fwrite(both_gene, output_path, append = FALSE, col.names = TRUE)
output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/ann_gene.csv"
fwrite(ann_gene, output_path, append = FALSE, col.names = TRUE)
output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/gordon_gene.csv"
fwrite(gordon_gene, output_path, append = FALSE, col.names = TRUE)

gordon_combined_genes <- rbind(combined_gene, gordon)
output_path <- "/Users/adams/Documents/master/master 2/Stage en thesis/Data/sars_cov_2/HCIP/gordon_combined_genes.csv"
fwrite(gordon_combined_genes, output_path, append = FALSE, col.names = TRUE)

gordon %>% pull(idBait) %>% unique
