
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


contrasts = readRDS("data/contrasts.hypertrophy.rds")

contrast_df = contrasts%>% mutate(contrast_id = paste(model, modal, tp, sep= "_"))%>%
  rename(gene= MgiSymbol)%>%
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
  mutate(contrast_id = paste("Hs_singlecell",Comparison, CellType, sep= "_"),
         FDR= 10^-log10.P.) %>% 
  rename(gene= Gene)%>%
  dplyr::select(contrast_id, gene, logFC, FDR)

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

### combine

joint_contrast_df= rbind(contrast_df, reheat.r, sc.gex2, contrast_m, PE)

joint_contrast_df <-
joint_contrast_df %>% group_by(contrast_id) %>% 
  mutate(ranks1= sign(logFC)*-log10(FDR),
         ranks2= rank(-ranks1),
         max.ranks= max(ranks2), 
         ranks3= ranks2/max.ranks, 
         gene = toupper(gene))


unique(joint_contrast_df$contrast_id)


saveRDS(joint_contrast_df, "data/contrasts_query_df.rds")
joint_contrast_df= readRDS( "data/contrasts_query_df.rds")

get_top_consistent_gene<-
  function(joint_contrast_df,
           query_contrasts= c("tac_ribo_2wk", "fetal_rna_fetal1", "HCMvsNF_Fibroblast"), 
           alpha= 0.05, 
           cutoff= 15){
    
    
    contrast_df_filt= joint_contrast_df %>% 
      filter(contrast_id %in% query_contrasts)%>%
      filter(FDR< alpha)
    
    venn= table(contrast_df_filt$gene)
    p.hist= enframe(venn) %>% ggplot(., aes(x= factor(value)))+
      geom_histogram(stat="count")+
      labs(x= "number of contrasts reporting gene")
    
    intersect_genes= names(venn[venn == length(query_contrasts)])
    
    # calculate the median normalized rank across contrasts
    df.median= joint_contrast_df %>%
      filter(contrast_id %in% query_contrasts,
                                gene %in% intersect_genes)%>% 
      group_by(gene)%>%
      summarise(m.r= median(ranks3))
    
    df.full= df.median %>%
      arrange(desc(m.r)) %>%
      left_join(contrast_df_filt%>%
                  select(gene, logFC, FDR, contrast_id), by= "gene")
    
    top_dn= df.median%>% arrange(desc(m.r))%>% slice(1:cutoff)%>% pull(gene)
    top_up= df.median%>% arrange(m.r)%>% slice(1:cutoff)%>% pull(gene)
    
    
     p.top_genes= 
       ggplot(df.full %>% filter(gene %in% c(top_dn, top_up))%>% 
             mutate(gene= factor(gene, levels= c(top_dn, top_up))), 
           aes(x= gene, y= logFC))+
      geom_boxplot(outlier.colour = NA)+
      geom_jitter(aes( color= contrast_id))+
      theme(axis.text.x = element_text(angle= 60 , hjust= 1))+
      geom_hline(yintercept = 0)
     
     return(list(p.hist=p.hist, 
                 genes= intersect_genes, 
                 df= df.full, 
                 p.top_genes= p.top_genes
                 ))
  }


res= get_top_consistent_gene(joint_contrast_df,
                             query_contrasts = c("DCMvsHCM_Fibroblast","swim_rna_2wk"),
                            cutoff= 10 )

res$p.hist
res$df
res$p.top_genes
