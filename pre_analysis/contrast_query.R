
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
  rbind(fetal1, fetal2)%>% 
  mutate( ## change here for adequate homologue mapping
         GeneDescription= "")%>%
  rename(
         PValue= P.Value,
         FDR= adj.P.Val)%>% 
  dplyr::select(logFC, PValue, FDR, gene, GeneDescription, tp,modal,model)

contrasts= rbind(fetal.df, mm2)

#contrasts = readRDS("data/contrasts.hypertrophy.rds")

contrast_df = contrasts%>% mutate(contrast_id = paste(model, modal, tp, sep= "_"))%>%
  dplyr::select(contrast_id, gene, logFC, FDR)%>%
  mutate(contrast_id= str_replace_all(contrast_id, "_fetal_rna", ""), 
         contrast_id = ifelse(!grepl("fetal", contrast_id), paste0("Mm_", contrast_id), contrast_id),
         contrast_id= str_replace_all(contrast_id, "fetal_rna_", ""))
  

unique(contrast_df$contrast_id)

#ReHeaT
ranks = readRDS("data/study_ranks.rds")
reheat.r= ranks%>% mutate(contrast_id= "Hs_ReHeaT")%>%
  rename(FDR= fisher_pvalue, 
         logFC= mean_lfc)%>% 
  dplyr::select(contrast_id, gene, logFC, FDR)

#chaffin 
sc.gex= read.csv("data/sc_gex_chaffin.csv")%>% 
  as_tibble()%>%
  mutate(Significant= factor(ifelse(Significant==1, TRUE, FALSE)))

sc.gex2= sc.gex %>% 
  mutate(contrast_id = paste("Hs_sc",Comparison, CellType, sep= "_"),
         FDR= 10^-log10.P., 
         logFC = ifelse(grepl("DCMvsHCM", contrast_id), logFC*-1, logFC))   %>% 
  rename(gene= Gene)%>%
  dplyr::select(contrast_id, gene, logFC, FDR)

#flip direction of contrast comp (logFC flipped above)
sc.gex2$contrast_id= str_replace_all(pattern = "DCMvsHCM",replacement =  "HCMvsDCM",sc.gex2$contrast_id)

## change the directionality of the hcm dcm comparison:
sc.gex2  = sc.gex2%>% 
  mutate(logFC = ifelse(grepl("DCMvsHCM", contrast_id), logFC * -1, logFC),
         contrast_id= str_replace_all(contrast_id, "DCMvsHCM", "HCMvsDCM"))

#magnet

magnet= readRDS("raw_data/magnet/contrast_list.rds")

contrast_m= lapply(names(magnet), function(x){
  magnet[[x]] %>% 
    mutate(contrast_id= paste0("Hs_bulk_", x))%>%
    rename(FDR= adj.P.Val)%>%
    dplyr::select(contrast_id, gene, logFC, FDR)
})%>% do.call(rbind, .)

#PE

PE= readRDS("raw_data/invitro_PE/contrast_df.rds")%>%
  dplyr:: select(contrast_id, gene, logFC, FDR)%>%
  mutate(gene =toupper(gene))


## proteome

swim_prot= read_csv("raw_data/proteome/Maus Swim/Limma_results_V2.csv")
swim_prot <- swim_prot %>% 
  mutate(gene = sapply(strsplit(gene_name, "\\|"), `[`, 1) %>% 
           sub("\\..*", "", .)) %>%
  mutate(contrast_id = "Mm_swim_prot", 
         contrast_id = ifelse(grepl("2d", sample), 
                              paste0(contrast_id,"_2d"),
                              paste0(contrast_id,"_2wk")))%>%
  rename(FDR = fdr.limma)%>%
  select(contrast_id,  gene,logFC, FDR)


tac_prot= read_csv("raw_data/proteome/Maus TAC/Limma_results_V1.csv")
tac_prot <- tac_prot %>% 
  mutate(gene = sapply(strsplit(gene_name, "\\|"), `[`, 1) %>% 
           sub("\\..*", "", .)) %>%
  mutate(contrast_id = "Mm_tac_prot", 
         contrast_id = ifelse(grepl("2d", sample), 
                              paste0(contrast_id,"_2d"),
                              paste0(contrast_id,"_2wk")))%>%
  rename(FDR = fdr.limma)%>%
  select(contrast_id,  gene,logFC, FDR)

hs_prot= read_csv("raw_data/proteome/Human Prosser Paper/Limma_results_V1.csv")

hs_prot <- hs_prot %>% 
  mutate(gene = sapply(strsplit(Gene.names, ";"), `[`, 1)) %>%
  mutate(contrast_id = comparison.label)%>%
  rename(FDR = fdr.limma)%>%
  select(contrast_id,  gene,logFC, FDR)
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


joint_contrast_df= rbind(contrast_df, 
                         reheat.r, 
                         sc.gex2, 
                         contrast_m, 
                         PE,
                         swim_prot,
                         tac_prot,
                         hs_prot)

joint_contrast_df <-
joint_contrast_df %>% 
  group_by(contrast_id) %>% 
  mutate(min.FDR = min(FDR[FDR > 0]),  ## modify FDR values that are INF to smallest non-zero FDR within that contrast
         FDR_mod= ifelse(FDR== 0, min.FDR, FDR), # replace 
         new_weight= sign(logFC)*-log10(FDR_mod))%>% ## calculate a weight based on signed -log10FDR 
  mutate(ranks1= sign(logFC)*-log10(FDR),
         ranks2= rank(-ranks1),
         max.ranks= max(ranks2), 
         ranks3= ranks2/max.ranks, 
         gene = toupper(gene), # this will be updated in species_translation.R
         sig= FDR<0.05, 
         sig= ifelse(grepl("Rn|Mm", contrast_id) & logFC<log2(1.5), FALSE, sig)) # as discussed we add this 


unique(joint_contrast_df$contrast_id)

joint_contrast_df = 
  joint_contrast_df %>% 
  mutate(cc= ifelse(grepl("Rn|Mm", contrast_id), "A", "C") ,
         cc= ifelse(grepl("fetal", contrast_id), "D", cc),
         cc= ifelse(grepl("HCM", contrast_id), "B", cc))

saveRDS(joint_contrast_df, "data/contrasts_query_df_untranslated.rds2.rds")



