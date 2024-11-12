
# translate symbols -------------------------------------------------------

#there are different conversions necessary between mouse, rat, human and genes and proteins.

#1. Mouse_protein -> Mouse_gene 
#2. Rat_protein -> Rat_gene
#2. Mouse_gene -> Human_gene
#3. Rat_gene -> Human_gene

# we will use biomaRt to perform each conversion, and we will quantify the
# unmapped as well as the double mapped features. 


# Load packages
library(biomaRt)
library(tidyverse)
library(cowplot)

project_path = "~/R-projects/Collaborations/cheerio/"

joint_contrast_df<-readRDS(paste0(project_path,"data/contrasts_query_df_untranslated2.rds" ))

unique(joint_contrast_df$contrast_id)

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
rat = useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")


#listDatasets(rat)%>% filter(grepl("r", dataset))
# 1. Mouse_protein -> Mouse_gene --------------------------------------

protein_ids <- joint_contrast_df%>% filter(grepl("ENSMUSP", gene))%>% pull(gene)%>% unique()

mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

# Get the mapping between protein IDs and gene symbols
mapped_genes <- biomaRt::getBM(
  attributes = c("ensembl_peptide_id","ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_peptide_id",
  values = protein_ids,
  mart = mouse
)

write.csv(mapped_genes, paste0(project_path, "data/translated_mm_prot_gene.csv"))


# 2. Mouse_gene -> human_gene ---------------------------------------------

mouse_gene_ids <- joint_contrast_df %>% filter(grepl("Mm", contrast_id), !grepl("ENSMUSP", gene)) %>% pull(gene) %>% unique()

# Map Mouse gene IDs to Human genes
mouse_to_human <- biomaRt::getBM(
  attributes = c("ensembl_gene_id","external_gene_name", "hsapiens_homolog_ensembl_gene", "hsapiens_homolog_associated_gene_name"),
  filters = "external_gene_name",
  values = mouse_gene_ids,
  mart = mouse
)
write.csv(mouse_to_human, paste0(project_path, "data/translated_mouse_human_gene.csv"))

# Quantify unmapped and duplicated
#mouse_human_stats <- quantify_mappings(mouse_to_human, mouse_gene_ids, "ensembl_gene_id")


# 3. Rat Gene -> Human Gene -----------------------------------------
rat_gene_ids <- joint_contrast_df %>% filter(grepl("Rn", contrast_id) & grepl("ENSRN", gene))%>% pull(gene) %>% unique()
rat_gene_ids_2 <- joint_contrast_df %>% filter(grepl("Rn", contrast_id) & !grepl("ENSRN", gene))%>% pull(gene) %>% unique()
# Map Rat gene IDs to Human genes
rat_to_human <- biomaRt::getBM(
  attributes = c("external_gene_name","ensembl_gene_id", "hsapiens_homolog_ensembl_gene", 
                 "hsapiens_homolog_associated_gene_name"),
  #filters = "external_gene_name", # doublecheck this
  filters = "ensembl_gene_id",
  values = rat_gene_ids,
  mart = rat
)

rat_to_human2 <- biomaRt::getBM(
  attributes = c("external_gene_name","ensembl_gene_id", "hsapiens_homolog_ensembl_gene", 
                 "hsapiens_homolog_associated_gene_name"),
  #filters = "external_gene_name", # doublecheck this
  filters = "ensembl_gene_id",
  values = rat_gene_ids_2,
  mart = rat
)

rat_to_human <- bind_rows(rat_to_human2, rat_to_human)
# Quantify unmapped and duplicated
#rat_human_stats <- quantify_mappings(rat_to_human, rat_gene_ids, "ensembl_gene_id")

write.csv(rat_to_human, paste0(project_path, "data/translated_rat_human_gene.csv"))


# assess mappings and update the joint contrast_df -------------------------------------
m_prot_to_m_gene <-read_csv(paste0(project_path, "data/translated_mm_prot_gene.csv"))
mouse_to_human<- read_csv(paste0(project_path, "data/translated_mouse_human_gene.csv"))
rat_to_human <- read_csv(paste0(project_path, "data/translated_rat_human_gene.csv"))

table(rat_gene_ids_2 %in% rat_to_human$external_gene_name)
table(rat_gene_ids %in% rat_to_human$ensembl_gene_id)

quantify_mappings <- function(mapped_df, original_ids, id_column, id_column_map) {
  mapped_df <- mapped_df %>% distinct(!!sym(id_column), !!sym(id_column_map))
  unmapped <- setdiff(original_ids, mapped_df[[id_column]])
  one_to_one <- mapped_df %>% group_by_at(id_column) %>% filter(n() ==1) 
  one_to_many <- mapped_df %>% group_by_at(id_column) %>% filter(n() > 1)
  many_to_one <- mapped_df %>% group_by_at(id_column_map) %>% filter(n() > 1)
  list(
    one_to_none = (unmapped),
    one_to_many = one_to_many %>% pull(id_column),  
    many_to_one = many_to_one  %>% pull(id_column),
    one_to_one = one_to_one %>% pull(id_column)
  )
 
}

summarize_mapping_stats <- function(mapping_list) {
  lapply(mapping_list, length)
}

mapping_results <- list(
  mouse_protein_to_mouse_gene = quantify_mappings(m_prot_to_m_gene, protein_ids, "ensembl_peptide_id", "external_gene_name"),
  mouse_to_human = quantify_mappings(mouse_to_human, mouse_gene_ids, "external_gene_name", "hsapiens_homolog_associated_gene_name"),
  rat_to_human = quantify_mappings(rat_to_human, rat_gene_ids, "ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),
  rat_to_human_2 = quantify_mappings(rat_to_human, rat_gene_ids_2, "external_gene_name", "hsapiens_homolog_associated_gene_name")
)


# Compute statistics for each mapping
mapping_stats <- lapply(mapping_results, summarize_mapping_stats)

# Create a summary table for all mappings
create_summary_table <- function(stats) {
  tibble(
    Conversion = c("Mouse Protein to Mouse Gene", "Mouse Gene to Human Gene", "Rat Gene to Human Gene", "Rat Gene to Human Gene (External)"),
    one_to_none = sapply(stats, function(x) x$one_to_none),
    one_to_one = sapply(stats, function(x) x$one_to_one),
    one_to_many = sapply(stats, function(x) x$one_to_many),
    many_to_one = sapply(stats, function(x) x$many_to_one)
  )
}

mapping_summary <- create_summary_table(mapping_stats)

# Transform the summary table for easier plotting
mapping_summary_long <- mapping_summary %>%
  pivot_longer(cols = c(one_to_none, one_to_one, one_to_many, many_to_one), names_to = "Status", values_to = "Count") %>%
  mutate(logCount = ifelse(is.infinite(log10(Count)), 0, log10(Count))) %>%
  group_by(Conversion) %>%
  mutate(prop = Count / sum(Count))


palette <- rev(c("#2E86AB", "#72B5B7", "#F29E4C", "#A6339E"))

p.mapping.efficiency<- mapping_summary_long%>%
  mutate(Status = factor(Status, levels=rev(c("one_to_one", 
                                          "many_to_one", 
                                          "one_to_many", 
                                          "one_to_none"))))%>%
ggplot(., aes(x = Conversion, y =prop, fill = Status)) +
  geom_col(color="black" )+
  geom_text(aes(label = ifelse(prop >= 0.05, scales::percent(prop, accuracy = 1), "")), 
            position = position_stack(vjust = 0.5), # Center the labels within the bars
            color = "white") +
  scale_fill_manual(values = palette)+
  theme_cowplot()+
  scale_y_continuous(labels = scales::percent)+
  theme(axis.text.x = element_text(angle= 90, hjust= 1, vjust= 0.5), 
        axis.line = element_blank(), 
        axis.ticks.x = element_blank())+
  labs(x="", y="", fill ="")
p.mapping.efficiency

pdf("figures/mapping.efficiency.pdf",
    height= 5, width = 4)
p.mapping.efficiency
dev.off()

# update the contrast_df --------------------------------------------------
mapping_results$mouse_protein_to_mouse_gene
mouse_prot_mapped<- joint_contrast_df%>% 
  filter(grepl("prot", contrast_id))%>%
  left_join(m_prot_to_m_gene%>%
              filter(ensembl_peptide_id %in% mapping_results$mouse_protein_to_mouse_gene$one_to_one)%>%
              distinct(ensembl_peptide_id, external_gene_name),
            by= c("gene"= "ensembl_peptide_id"))%>% 
  drop_na()
mapping_results$rat_to_human
#map the rat
rat_mapped= joint_contrast_df%>% 
  filter(grepl("Rn", contrast_id) & !grepl("prot", contrast_id))%>%
  left_join(rat_to_human %>% 
              filter(ensembl_gene_id %in% mapping_results$rat_to_human$one_to_one)%>% 
              distinct(ensembl_gene_id, hsapiens_homolog_associated_gene_name),
            by= c("gene"= "ensembl_gene_id"))%>% 
  rename(gene_orig = gene, 
         gene = hsapiens_homolog_associated_gene_name)%>% 
  dplyr::select(contrast_id, gene, everything())%>% 
  drop_na()

rat_mapped2= joint_contrast_df%>% 
  filter(grepl("Rn", contrast_id) & grepl("prot", contrast_id))%>%
  left_join(rat_to_human %>% 
              filter(external_gene_name %in% mapping_results$rat_to_human_2$one_to_one)%>% 
              distinct(external_gene_name, hsapiens_homolog_associated_gene_name),
            by= c("gene"= "external_gene_name"))%>% 
  rename(gene_orig = gene, 
         gene = hsapiens_homolog_associated_gene_name)%>% 
  dplyr::select(contrast_id, gene, everything())%>% 
  drop_na()

#map the mouse
mouse_mapped<- joint_contrast_df%>% 
  filter(grepl("Mm", contrast_id))%>%
  left_join(mouse_to_human %>% 
              filter(external_gene_name %in% (mouse_to_human_l$one_to_one))%>% 
              distinct(external_gene_name, hsapiens_homolog_associated_gene_name),
            by= c("gene"= "external_gene_name"))%>% 
  rename(gene_orig = gene, 
         gene = hsapiens_homolog_associated_gene_name)%>% 
  dplyr::select(contrast_id, gene, everything())%>% 
  drop_na()

mouse_prot_mapped_to_gene<-mouse_prot_mapped%>% 
  left_join(mouse_to_human %>% 
              filter(external_gene_name %in% unlist(mouse_to_human_l$one_to_one))%>% 
              distinct(external_gene_name, hsapiens_homolog_associated_gene_name),
            by= c("external_gene_name"= "external_gene_name"))%>% 
  dplyr::select(-gene)%>%
  rename(gene_orig = external_gene_name, 
         gene = hsapiens_homolog_associated_gene_name)%>% 
  dplyr::select(contrast_id, gene,gene_orig,  everything() )%>% 
  drop_na()

Animal_df <- rbind(rat_mapped,rat_mapped2, mouse_mapped, mouse_prot_mapped_to_gene )

mapped_contrasts<- unique(Animal_df$contrast_id )
#make sure all animal contrasts have been mapped->
unique(joint_contrast_df$contrast_id)[!unique(joint_contrast_df$contrast_id) %in% mapped_contrasts]

Human_df<- joint_contrast_df %>%
  filter(grepl("Hs", contrast_id))%>%
  mutate(gene_orig= gene)%>%
  dplyr::select(contrast_id, gene,gene_orig,  everything() )
      
joint_df_translated<- rbind(Animal_df, Human_df)   
unique(joint_df_translated$contrast_id)

joint_df_translated%>% 
  group_by(contrast_id)%>%
  summarize(unique_gene_count = n_distinct(gene))%>%
  pull(unique_gene_count)%>% hist(., breaks=100)
  

saveRDS(joint_df_translated, "data/contrasts_query_df_translated2.rds")
