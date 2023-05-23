contrasts = readRDS("data/contrasts.hypertrophy.rds")

sc.gex= read.csv("data/sc_gex_chaffin.csv")%>% 
  as_tibble()%>%
  mutate(Significant= factor(ifelse(Significant==1, TRUE, FALSE)))


directed_signature = readRDS("data/signature.rds")

directed_signature

contrast_df = contrasts%>% mutate(contrast_id = paste(model, modal, tp, sep= "_"))%>%
  rename(gene= MgiSymbol)%>%
  select(contrast_id, gene, logFC, FDR)

reheat.r= ranks%>% mutate(contrast_id= "ReHeaT")%>%
  rename(FDR= fisher_pvalue, 
         logFC= mean_lfc)%>% 
  select(contrast_id, gene, logFC, FDR)

sc.gex2= sc.gex %>% 
  mutate(contrast_id = paste(Comparison, CellType, sep= "_"),
         FDR= 10^-log10.P.) %>% 
  rename(gene= Gene)%>%
  select(contrast_id, gene, logFC, FDR)

joint_contrast_df= rbind(contrast_df, reheat.r, sc.gex2)

joint_contrast_df <-
joint_contrast_df %>% group_by(contrast_id) %>% 
  mutate(ranks1= sign(logFC)*-log10(FDR),
         ranks2= rank(-ranks1),
         max.ranks= max(ranks2), 
         ranks3= ranks2/max.ranks, 
         gene = toupper(gene))


unique(joint_contrast_df$contrast_id)

saveRDS(joint_contrast_df, "data/contrasts_query_df.rds")

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
