library(data.table)
library(tidyverse)
library(taxize)
library(myTAI)

# Load data -------------------------------------------------------------------------------------------------------------------------

#	Load protein-protein interactions
tbl_VirHostNet_reanalysis <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/Virhostnet/VirHostNet_reanalysis.tsv") %>% as_tibble

#	Load VirHostNet file
tbl_VirHostNet_taxid <- fread("/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/Virhostnet/VirHostNet_taxid.tsv") %>% as_tibble

#	Create dataframe with human proteins interacting with viral proteins
tbl_Vir1 <- tbl_VirHostNet_taxid %>% 
	filter(str_detect(Genename1, "HUMAN")) %>%
	separate(UniprotID1, into = c("uni1", "humanID")) %>%
	separate(UniprotID2, into = c("uni2", "viralID")) %>%
	separate(Genename1, into = c("humangene", "humanspecies"), sep = "_") %>% 
	separate(Genename2, into = c("viralgene", "viralspecies"), sep = "_") %>% 
	mutate(humantaxid = taxid1) %>%
	mutate(viraltaxid = taxid2) %>%
	select(ends_with("ID") | ends_with("gene") | ends_with("taxid") | ends_with("species"))

tbl_Vir2 <- tbl_VirHostNet_taxid %>% 
	filter(str_detect(Genename2, "HUMAN")) %>%
	separate(UniprotID2, into = c("uni1", "humanID")) %>%
	separate(UniprotID1, into = c("uni2", "viralID")) %>%
	separate(Genename2, into = c("humangene", "humanspecies"), sep = "_") %>% 
	separate(Genename1, into = c("viralgene", "viralspecies"), sep = "_") %>% 
	mutate(humantaxid = taxid2) %>%
	mutate(viraltaxid = taxid1) %>%
	select(ends_with("ID") | ends_with("gene") | ends_with("taxid") | ends_with("species"))

tbl_VirHostNet_human <- rbind(tbl_Vir1, tbl_Vir2) %>%
	filter(!viraltaxid == humantaxid)

#	Map virus name to taxid
df_virhostnet_analysis <- merge(tbl_VirHostNet_human, tbl_VirHostNet_reanalysis, by.x = "humanID", by.y = "Prey") %>% 
	as_tibble() %>%
	separate(viraltaxid, into = c("tax", "taxid"), sep = ":", remove = FALSE)

taxids <- df_virhostnet_analysis %>% separate(viraltaxid, into = c("tax", "taxid"), sep = ":") %>% pull(taxid)
taxon_summary <- ncbi_get_taxon_summary(id = taxids) %>% as_tibble

#	Map interacting proteins to viral proteins
df_virhostnet_analysis_names <- df_virhostnet_analysis %>% 
	merge(taxon_summary, by.x = "taxid", by.y = "uid") %>%
	unique()

# Get numbers -------------------------------------------------------------------------------------------------------------------------

#	375 host proteins interact with SARS-CoV-2
tbl_VirHostNet_reanalysis %>% 
	select(Prey) %>% unique()

#	328 host proteins in other viruses
df_virhostnet_analysis_names %>% as_tibble %>% 
	filter(!str_detect(name, "Severe acute respiratory")) %>% 			# Coronavirus?
	select(humanID) %>% unique()

#	Load genes background
fasta_headers <- "/Users/adams/Documents/PhD/SARS-CoV-2/Data/Workspace/fasta/fasta_headers.csv"
unmapped_headers <- fread(fasta_headers, header = FALSE) %>% as_tibble
background <- unmapped_headers %>% separate(V3, into=c("pre", "gene"), sep = (" GN=")) %>% drop_na(gene) %>%
	separate(gene, into = c("gene_name", "post", "ppost") ,sep=(" ")) %>%
	rename(humanID = V2) %>%
	unique()

#	69320 proteins that don't interact with SARS-CoV-2
tbl_background <- background %>% 
	filter(!humanID %in% tbl_VirHostNet_reanalysis$Prey) 

tbl_background %>%
	select(humanID) %>%
	unique()

#	7089 proteins in VirHostNet (excluding SARS1 and SARS2)
merge(tbl_VirHostNet_human, tbl_background, by = "humanID") %>% 
	as_tibble() %>%
	separate(viraltaxid, into = c("tax", "taxid"), sep = ":", remove = FALSE) %>%
	merge(taxon_summary, by.x = "taxid", by.y = "uid") %>%
	filter(!str_detect(name, "Severe acute respiratory")) %>% 			# Coronavirus?
	select(humanID) %>%
	as_tibble() %>%
	unique()

#	Fisher test ---------------------------------------------------------------------------------------------------------------
m <- 7185	# In VirHostNet
n <- 62510	# Not in VirHostNet
k <- 375	# Interact with SARS
x <- 44	# Both in VirHostNet and SARS

probability <- dhyper(x, m, n, k, log = FALSE)

m <- 375	# In VirHostNet
n <- 62510	# Not in VirHostNet
k <- 375	# Interact with SARS
x <- c(0:375)	# Both in VirHostNet and SARS

probability <- dhyper(x, m, n, k, log = FALSE)