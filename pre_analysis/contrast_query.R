
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
  select(contrast_id,  gene,logFC, pval, FDR)


tac_prot= read_csv("raw_data/proteome/Maus TAC/Limma_results_V1.csv")
tac_prot <- tac_prot %>% 
  mutate(gene = sapply(strsplit(gene_name, "\\|"), `[`, 1) %>% 
           sub("\\..*", "", .)) %>%
  mutate(contrast_id = "Mm_tac_prot", 
         contrast_id = ifelse(grepl("2d", sample), 
                              paste0(contrast_id,"_2d"),
                              paste0(contrast_id,"_2wk")))%>%
  rename(FDR = fdr.limma,
         pval = pvalue.limma)%>%
  select(contrast_id,  gene,logFC, pval, FDR)

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
                         sc.gex2, 
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

unique(joint_contrast_df$contrast_id)

## add contrast category
joint_contrast_df = 
  joint_contrast_df %>% 
  mutate(cc= ifelse(grepl("Rn|Mm", contrast_id), "A", "C") ,
         cc= ifelse(grepl("fetal", contrast_id), "D", cc),
         cc= ifelse(grepl("HCM", contrast_id), "B", cc))

saveRDS(joint_contrast_df, "data/contrasts_query_df_untranslated2.rds")




# create function for fisher p test ---------------------------------------
source("sub/global.R")
source("sub/helper.R")

get_top_consistent_gene<-
  function(joint_contrast_df,
           query_contrasts= c("Mm_tac_ribo_2wk", "Hs_bulk_HCMvsNF"), 
           alpha= 0.05, 
           cutoff= 15,
           missing_prop= 10){
    
    
    contrast_df_filt= joint_contrast_df %>% 
      filter(contrast_id %in% query_contrasts)%>%
      filter(FDR< alpha)
    
    # get overview of intersect
    venns= table(contrast_df_filt$gene)
    
    x= split(contrast_df_filt$gene, contrast_df_filt$contrast_id)
    x= lapply(x, unique)
    
    p.venn= plot(euler(x, shape = "ellipse"), quantities = TRUE)
    
    intersect_genes= names(venns[venns >= (length(query_contrasts)*missing_prop/100)])
    
    df.msign= contrast_df_filt %>%
      dplyr::select(gene, contrast_id, logFC)%>% 
      filter(gene %in% intersect_genes) %>%
      group_by(gene)%>% 
      summarise(m.sign = mean(sign(logFC)))%>%
      mutate(top_ = ifelse(m.sign== 1, "upregulated", 
                           ifelse(m.sign ==-1 , "downregulated", "inconsistent")))
    
    top_up = df.msign %>%
      filter(top_== "upregulated")%>% pull(gene)
    
    top_dn = df.msign %>%
      filter(top_== "downregulated")%>% pull(gene)
    
    p.bar.intersect=
      ggplot(df.msign, aes(x= factor(top_, levels = c("upregulated", 
                                                      "downregulated", 
                                                      "inconsistent"))))+
      geom_bar()+
      labs(x= "consistency of regulation")+
      theme(axis.text.x = element_text(angle= 60, hjust= 1))
    
    p.int=  cowplot::plot_grid( p.venn,p.bar.intersect,
                                ncol = 2, labels = c("A", 
                                                     "B"), 
                                rel_widths = c(1,0.3))
    
    # calculate the median normalized rank across contrasts
    df.median= joint_contrast_df %>%
      filter(contrast_id %in% query_contrasts,
             gene %in% intersect_genes)%>% 
      group_by(gene)%>%
      summarise(m.r= mean(ranks3))
    
    df.full= df.median %>%
      arrange(desc(m.r)) %>%
      left_join(contrast_df_filt%>%
                  dplyr::select(gene, logFC, FDR, contrast_id), by= "gene")
    
    if(length(intersect_genes)!= 0){
      top_dn2= df.median%>% filter(gene %in% top_dn) %>% arrange(desc(m.r))%>% slice(1:cutoff)%>% pull(gene)
      top_up2= df.median%>% filter(gene %in% top_up) %>% arrange(m.r)%>% slice(1:cutoff)%>% pull(gene)
      
      p.top_gene_up=
        df.full %>% 
        filter(gene %in% c( top_up2))%>% 
        mutate(gene= factor(gene, levels= c(top_up2)))%>%
        ggplot(., 
               aes(x= gene, y= logFC))+
        geom_boxplot(outlier.colour = NA)+
        geom_jitter(aes( color= contrast_id))+
        theme(axis.text.x = element_text(angle= 60 , hjust= 1))
      #geom_hline(yintercept = 0)
      
      p.top_gene_dn=
        df.full %>% 
        filter(gene %in% c(top_dn2 ))%>% 
        mutate(gene= factor(gene, levels= c(top_dn2 )))%>%
        ggplot(., 
               aes(x= gene, y= logFC))+
        geom_boxplot(outlier.colour = NA)+
        geom_jitter(aes( color= contrast_id))+
        theme(axis.text.x = element_text(angle= 60 , hjust= 1))
      #geom_hline(yintercept = 0)      
      
      p.top_genes = cowplot::plot_grid(p.top_gene_up, p.top_gene_dn, 
                                       ncol = 1,
                                       labels= c("top upregulated", 
                                                 "top downregulated"))
      
      ##add hmap
      plot.genes= unique(c(top_dn, top_up))
      
      mat= joint_contrast_df %>% 
        dplyr::select(contrast_id, gene,  logFC)%>%
        filter(gene %in% plot.genes,
               contrast_id %in% query_contrasts)%>%
        pivot_wider(names_from = contrast_id, values_from = logFC, values_fn = mean)%>%
        as.data.frame()%>%
        filter(!is.na(gene))%>%
        column_to_rownames("gene")
      col_names_plot= c(top_dn2, top_up2)
      
      na_sums_per_row <- apply(mat, 1, function(row) sum(is.na(row)))
      mat_subset= mat[na_sums_per_row < 0.5 *  ncol(mat),]
      
      x= rownames(mat_subset)
      x[!x %in% col_names_plot] <- ""
      
      hmap_top <- ComplexHeatmap::Heatmap(t(mat_subset),
                                          #rect_gp = gpar(fill = "grey", lwd = 1),
                                          name = "logFC", 
                                          na_col = "black",
                                          border_gp = gpar(col = "black", lty = 1),
                                          cluster_columns = T,
                                          cluster_rows= T, 
                                          
                                          column_labels = x,
                                          column_names_gp = gpar(fontsize = 9),
                                          row_names_side = "left",
                                          row_dend_side = "left"
                                          
      )
      
    }else{ 
      p.top_genes= NULL
      hmap_top <- NULL
    }
    
    ## add hmap= 
    
    
    
    
    return(list(p.hist=p.int, 
                genes= list("i"= intersect_genes, 
                            "u"= top_up, 
                            "d"= top_dn),
                df= df.full, 
                p.top_genes= p.top_genes,
                hmap_top= hmap_top
    ))
  }


run_fisher_meta = function(meta_list, n_missing = 3){
  library(survcomp)
  # Getting p-values from limma
  limma_pvals = get_all_limma(meta_list = meta_list, "adj.P.Val")
  
  # Use only genes that are present in all experiments (missing in n at most)
  limma_results_mat = limma_pvals[rowSums(is.na(limma_pvals))<=n_missing,]
  
  # Fisher combined test
  fisher_pvals = apply(limma_results_mat, 1, function(x){ 
    survcomp::combine.test(x, "fisher", na.rm = T)
  })
  
  fisher_pvals_adj = sort(p.adjust(fisher_pvals,"BH"))
  
  return(fisher_pvals_adj)
}


df.full %>%
  select(gene, FDR)
