h# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2023 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2023-12-18
#
# Script Name:    ~/R-projects/Collaborations/shiny_hypertophy/species_translation.R
#
# Script Description:
# translate between species.  

library(msigdbr)
library(annotables)

source("utils_pre_analysis.R")
contrast_df = readRDS("data/contrasts_query_df.rds")


rat_genes= contrast_df %>%
  filter(grepl("Rn", contrast_id))%>% pull(gene) %>% unique()%>% str_to_title()

mouse_genes= contrast_df %>%
  filter(grepl("Mm", contrast_id))%>% pull(gene) %>% unique()%>% str_to_title()


rat_hs <- translate_species_to_hs("Rattus norvegicus",rat_genes)
mm_hs <- translate_species_to_hs("Mus musculus",mouse_genes)


# update the contrast df 
contrast_df3 <-contrast_df %>% 
  mutate( gene_orig= ifelse(!grepl("Hs", contrast_id), str_to_title(gene), gene))

contrasts_mouse= contrast_df3%>%
  filter(grepl("Mm", contrast_id))%>%
  filter(gene_orig %in% mm_hs$genes$oto)%>%
  dplyr::select(-gene)%>%
  left_join(mm_hs$df%>%
              distinct(gene_symbol, human_gene_symbol)%>%
              rename(gene_orig= gene_symbol,
                     gene= human_gene_symbol))
  

contrasts_rat= contrast_df3%>%
  filter(grepl("Rn", contrast_id))%>%
  filter(gene_orig %in% rat_hs$genes$oto)%>%
  dplyr::select(-gene)%>%
  left_join(rat_hs$df%>%
              distinct(gene_symbol, human_gene_symbol)%>%
              rename(gene_orig= gene_symbol,
                     gene= human_gene_symbol))

contrasts_h= contrast_df %>% 
  filter(grepl("Hs", contrast_id))%>% 
  mutate(gene_orig= gene)
  
contrasts_df_updated= rbind(rbind(contrasts_h,contrasts_rat[,colnames(contrasts_h)]),
      contrasts_mouse[,colnames(contrasts_h)]
      )

saveRDS(contrasts_df_updated, "app_data/contrasts_query_df.rds")

