
## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-05-26
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## collect all contrasts and provide them together in tidy format

library(tidyverse)

## A. MM models
mm <- readRDS("data/contrasts_mm_limma.rds")
mm2= mm %>% mutate(GeneDescription= "",
              #MgiSymbol = gene, 
              FDR= adj.P.Val,
              PValue= P.Value)%>%
  dplyr::select(logFC, PValue, FDR, gene, GeneDescription, tp,modal,model
                )

## add fetal expression data:
# fetal1= readRDS(file = "data/fetalDEgenes_GSE52601.rds") %>% 
#   mutate(tp= "Hs_fetal_Akat14", 
#          modal= "rna",
#          model ="fetal") 

fetal2= readRDS(file = "data/fetalDEgenes_PRJNA522417.rds")%>% 
  mutate(tp= "Hs_fetal_Spurell22", 
         modal= "rna",
         model ="fetal") 

fetal.df= 
  #rbind(fetal1, fetal2)%>% 
  fetal2%>%
  mutate( ## change here for adequate homologue mapping
         GeneDescription= "")%>%
  rename(
         PValue= P.Value,
         FDR= adj.P.Val)%>% 
  dplyr::select(logFC, PValue, FDR, gene, GeneDescription, tp,modal,model)

contrasts= rbind(fetal.df, mm2)

#contrasts = readRDS("data/contrasts.hypertrophy.rds")

contrast_df = contrasts%>% mutate(contrast_id = paste(model, modal, tp, sep= "_"))%>%
  dplyr::select(contrast_id, gene, logFC, PValue, FDR)%>%
  mutate(contrast_id= str_replace_all(contrast_id, "_fetal_rna", ""), 
         contrast_id = ifelse(!grepl("fetal", contrast_id), paste0("Mm_", contrast_id), contrast_id),
         contrast_id= str_replace_all(contrast_id, "fetal_rna_", ""))%>%
  rename(pval= PValue)
  

unique(contrast_df$contrast_id)

#ReHeaT
ranks = readRDS("data/study_ranks.rds")
reheat.r= ranks%>% mutate(contrast_id= "Hs_ReHeaT")%>%
  rename(FDR= fisher_pvalue, 
         logFC= mean_lfc)%>% 
  dplyr::select(contrast_id, gene, logFC, FDR)

#chaffin 
sc.gex <- read_csv("data/chaffinetal_degs.csv")%>%
  rename(gene= Gene)%>%
  mutate(contrast_id = paste("Hs_sc",Comparison, Cell_type, sep= "_"))%>%
  dplyr::select(contrast_id, gene, logFC, pval , FDR)
# 
# sc.gex= read.csv("data/sc_gex_chaffin.csv")%>% 
#   as_tibble()%>%
#   mutate(Significant= factor(ifelse(Significant==1, TRUE, FALSE)))

#magnet

magnet= readRDS("raw_data/magnet/contrast_list.rds")

contrast_m= lapply(names(magnet), function(x){
  magnet[[x]] %>% 
    mutate(contrast_id= paste0("Hs_bulk_", x))%>%
    rename(FDR= adj.P.Val, 
           pval = P.Value)%>%
    dplyr::select(contrast_id, gene, logFC, pval, FDR)
})%>% do.call(rbind, .)

# RAT in vitro 

x1= read_csv("../cheerio/raw_data/invitro_PE/Ribo_DEG_CPM10_NRVM_online.csv")
x2= read_csv("../cheerio/raw_data/invitro_PE/RNA_DEG_CPM10_NRVM_online.csv")

x1 = x1 %>% dplyr::rename(gene= Ensemble_ID)%>%
  dplyr::select(gene, logFC,PValue,  FDR)%>%
  rename(pval = PValue)%>%
  mutate(#contrast_field= "Rn_invitro", 
    contrast_id = "Rn_invitro_ribo")

x2 = x2 %>%
  dplyr::rename(gene= Ensemble_ID_RNA)%>%
  dplyr::select(gene, logFC,PValue,  FDR)%>%
  rename(pval = PValue)%>%
  mutate(#contrast_field= "Rn_invitro", 
    contrast_id = "Rn_invitro_rna")

rat_PE <- rbind(x1, x2)

## proteome
rat_prot= read_csv("raw_data/proteome/Rat_invitro/Limma_results_V1 PE vs Ctrl in vitro Kardios.csv" ,
                   skip = 1)
rat_prot<- rat_prot %>%
  rename(FDR = fdr.limma, 
         pval= pvalue.limma, 
         gene = Gene)%>%
  mutate(contrast_id= "Rn_invitro_prot")%>%
  select(contrast_id,  gene,logFC, pval, FDR)
  
swim_prot= read_csv("raw_data/proteome/Maus Swim/Limma_results_V2.csv")
swim_prot <- swim_prot %>% 
  mutate(gene = sapply(strsplit(gene_name, "\\|"), `[`, 1) %>% 
           sub("\\..*", "", .)) %>%
  mutate(contrast_id = "Mm_swim_prot", 
         contrast_id = ifelse(grepl("2d", sample), 
                              paste0(contrast_id,"_2d"),
                              paste0(contrast_id,"_2wk")))%>%
  rename(FDR = fdr.limma,
         pval = pvalue.limma)%>%
  dplyr::select(contrast_id,  gene,logFC, pval, FDR)

swim_prot
tac_prot= read_csv("raw_data/proteome/Maus TAC/Limma_results_V1.csv")
tac_prot <- tac_prot %>% 
  filter(comparison== "TAC - Sham_TAC")%>%
  mutate(gene = sapply(strsplit(gene_name, "\\|"), `[`, 1) %>% 
           sub("\\..*", "", .)) %>%
  mutate(contrast_id = "Mm_tac_prot", 
         contrast_id = ifelse(grepl("2d", sample), 
                              paste0(contrast_id,"_2d"),
                              paste0(contrast_id,"_2wk")))%>%
  rename(FDR = fdr.limma,
         pval = pvalue.limma)%>%
  dplyr::select(contrast_id,  gene,logFC, pval, FDR)
tac_prot%>% filter(gene=="ENSMUSP00000019683")
hs_prot= read_csv("raw_data/proteome/Human Prosser Paper/Limma_results_V1.csv")

hs_prot <- hs_prot %>% 
  mutate(gene = sapply(strsplit(Gene.names, ";"), `[`, 1)) %>%
  mutate(contrast_id = comparison.label)%>%
  rename(FDR = fdr.limma,
         pval = pvalue.limma)%>%
  select(contrast_id,  gene,logFC, pval, FDR)
hs_prot$contrast_id<-  str_replace_all(hs_prot$contrast_id,
                 pattern = "cHyp vs normal",
                 replacement = "Hs_bulk_prot_cHypvsNF" )
hs_prot$contrast_id<-  str_replace_all(hs_prot$contrast_id,
                                       pattern = "HCMrEF vs normal",
                                       replacement = "Hs_bulk_prot_HCMrEFvsNF" )

hs_prot$contrast_id<-  str_replace_all(hs_prot$contrast_id,
                                       pattern = "HCMpEF vs normal",
                                       replacement = "Hs_bulk_prot_HCMpEFvsNF" )
hs_prot$contrast_id<-  str_replace_all(hs_prot$contrast_id,
                                       pattern = "DCM vs normal",
                                       replacement = "Hs_bulk_prot_DCMvsNF" )
hs_prot<-hs_prot %>% filter(contrast_id %in% c("Hs_bulk_prot_cHypvsNF",
                                      "Hs_bulk_prot_HCMrEFvsNF",
                                      "Hs_bulk_prot_HCMpEFvsNF",
                                      "Hs_bulk_prot_DCMvsNF"))

unique(hs_prot$contrast_id)



# bind all data frames  ---------------------------------------------------

## join
joint_contrast_df= bind_rows(contrast_df, 
                         reheat.r, 
                         sc.gex, 
                         contrast_m, 
                         rat_PE,
                         rat_prot,
                         swim_prot,
                         tac_prot,
                         hs_prot)

## transform some o add rankings per contrasts
joint_contrast_df <-
joint_contrast_df %>% 
  group_by(contrast_id) %>% 
  mutate(min_FDR = min(FDR[FDR > 0]),  ## modify FDR values that are INF to smallest non-zero FDR within that contrast
         FDR_mod= ifelse(FDR== 0, min_FDR, FDR), # replace 0 values with the minimum
         signed_FDR_weight= sign(logFC)*-log10(FDR_mod))%>% ## calculate a weight based on signed -log10FDR 
  mutate(gene_rank= rank(-signed_FDR_weight),# a ranking on the signed pvalue
         n_genes= max(gene_rank), # corresponds to the number of genes
         gene_rank_norm= gene_rank/n_genes, #  normalizes the rank to number of genes
         sig= FDR<0.05, 
         sig= ifelse(grepl("Rn|Mm", contrast_id) & logFC<log2(1.5), FALSE, sig)) # we added this as an additional requirement for animal models to be significant

sort(unique(joint_contrast_df$contrast_id))

## add contrast category
joint_contrast_df = 
  joint_contrast_df %>% 
  mutate(cc= ifelse(grepl("Rn|Mm", contrast_id), "A", "C") ,
         cc= ifelse(grepl("fetal", contrast_id), "D", cc),
         cc= ifelse(grepl("HCM", contrast_id), "B", cc))


saveRDS(joint_contrast_df, "data/contrasts_query_df_untranslated2.rds")


