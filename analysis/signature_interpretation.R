# HEADER --------------------------------------------
#
# Author:     Jan D. Lanzer
# Copyright   Copyright 2024 - Jan D. Lanzer
# Email:      lanzerjan@gmail.com
#
# Date:     2024-11-15
#
# Script Name:    ~/R-projects/Collaborations/cheerio/analysis/signature_interpretation.R
#
# Script Description:
# interpret signatures

library(decoupleR)
library(tidyverse)
library(msigdbr)
cm_sets<- read_csv("data/prior_knowledge_cm.csv")
cm_sets


loadings.pca <- readRDS( "data/pcaloadings.rds")

pca.mat<-loadings.pca[,c("PC1", "PC2", 
                         #"PC3", "PC4",
                         "gene")]%>%
  as.data.frame()%>% 
  column_to_rownames("gene")%>%
  as.matrix()

hm<- msigdbr(species = "Homo sapiens", category = "H")%>%
  distinct(gs_name, gene_symbol)

C2<- msigdbr(species = "Homo sapiens", category = "C2")%>%
  distinct(gs_name, gene_symbol)

msig_res= decoupleR::run_ulm(mat = pca.mat, net= hm%>% distinct(),.source = gs_name, .target = gene_symbol 
                             )
msig_res= decoupleR::run_ulm(mat = pca.mat, net= C2,.source = gs_name, .target = gene_symbol 
)

write.csv(msig_res, "data/misg_enriched_in_pca_loadings.csv")

pways= msig_res %>%filter(statistic=="ulm") %>% mutate(p_adj= p.adjust(p_value))%>%
  group_by(condition, source)%>%
  filter(any(p_adj<0.01))%>%
  pull(source
  )

px= msig_res %>% mutate(p_adj= p.adjust(p_value)) %>%
  filter(source %in% pways)%>%
  ggplot(aes(x=score, y= reorder(source,score), fill = condition))+
  facet_grid(~ condition)+
  geom_col(aes(color= p_adj<0.01))+
  scale_color_manual(values= c("TRUE"= "black", "FALSE"= "white"))+
  theme_cowplot()+
  theme(axis.text.y = element_text(size=10))+
  geom_vline(xintercept = 0)
px
pdf("plots/big_pc_hmap.pdf_2",width= 20, height= 15)
px
dev.off()


msig_res_cm= decoupleR::run_ulm(mat = pca.mat, 
                                net= cm_sets,
                                .source = source, .target = target)


pways= msig_res_cm %>%filter(statistic=="ulm") %>% mutate(p_adj= p.adjust(p_value))%>%
  group_by(condition, source)%>%
  filter(any(p_adj<0.4))%>%
  pull(source
  )

px= msig_res_cm %>% mutate(p_adj= p.adjust(p_value)) %>%
  filter(source %in% pways)%>%
  ggplot(aes(x=score, y= reorder(source,score), fill = condition))+
  facet_grid(~ condition)+
  geom_col(aes(color= p_adj<0.01))+
  scale_color_manual(values= c("TRUE"= "black", "FALSE"= "white"))+
  theme_cowplot()+
  theme(axis.text.y = element_text(size=10))+
  geom_vline(xintercept = 0)
px
#plot also for a subset for main fig:
pways
pways_select  <- 
  c("NABA_CORE_MATRISOME", 
    "KEGG_FATTY_ACID_METABOLISM",                                                                                               
    "KEGG_OXIDATIVE_PHOSPHORYLATION" ,
    "KEGG_CITRATE_CYCLE_TCA_CYCLE",
    "HALLMARK_HYPOXIA" ,
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT" ,
    "REACTOME_MITOCHONDRIAL_TRANSLATION",
    "REACTOME_CROSSLINKING_OF_COLLAGEN_FIBRILS",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION",
    "HALLMARK_MTORC1_SIGNALING"
    
    
    
    
    
  )


px_select= msig_res %>% mutate(p_adj= p.adjust(p_value)) %>%
  filter(source %in% pways_select)%>%
  ggplot(aes(x=score, y= reorder(source,score), fill = condition))+
  facet_grid(~ condition)+
  geom_col(aes(color= p_adj<0.01))+
  scale_color_manual(values= c("TRUE"= "black", "FALSE"= "white"))+
  theme_cowplot()+
  theme(axis.text.y = element_text(size=10))+
  geom_vline(xintercept = 0)
px_select

pdf("plots/plot_msig_pca_loadings.pdf_2",
    width= 10, height= 3)
px_select
dev.off()

msig_res %>% mutate(p_adj= p.adjust(p_value)) %>%
  filter(source %in% pways, 
         condition== "PC2",
         score>0)%>% print(n=100)

#get human, rodent and joint signature:
df<- loadings.pca%>%as.data.frame() %>% rownames_to_column("gene")%>% as_tibble()
write.csv(loadings.pca, "data/pca_loadings_gene_list.csv")

h.genes <- df%>% 
  dplyr::filter(PC1<0.025)%>%
  arrange(desc(PC2))%>% 
  slice_head(n= 20)%>% 
  pull(gene)

m.genes <- df%>% 
  filter(PC1<0.025)%>%
  arrange((PC2))%>% 
  slice_head(n= 20)%>% 
  pull(gene)

joint<-df%>% 
  filter(PC1<0)%>%
  arrange(abs(PC2))%>% 
  slice_head(n= 20)%>% 
  pull(gene)


p.list= map( list(h.genes.sort, 
                  m.genes, 
                  joint), function(x){
                    gene<- x[1:8]
                    p1= rank.matrix[gene,] %>% 
                      as.data.frame()%>% rownames_to_column("gene") %>%
                      pivot_longer(-gene, names_to = "contrast_id", values_to= "logFC")%>%
                      left_join(c.df %>% 
                                  distinct(cc, contrast_id))%>%
                      filter(cc %in% c("A", "B"))%>%
                      ggplot(., aes(x= cc, y= logFC))+
                      facet_grid(~gene)+
                      geom_boxplot()+
                      geom_jitter(aes(color= contrast_id))+
                      geom_hline(yintercept = 0)+ theme_cowplot()
                  })

h.genes.sort <-   rank.matrix[h.genes,] %>%
  as.data.frame() %>% 
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "contrast_id", values_to= "logFC")%>%
  left_join(c.df %>% 
              distinct(cc, contrast_id))%>%
  filter(cc %in% c("A", "B"))%>%
  group_by(cc)%>%
  mutate(s.cc= sign(logFC))%>%
  group_by(gene, cc)%>%
  mutate(ss= sum(s.cc) )%>%
  distinct(gene, cc, ss)%>% 
  group_by(gene)%>%
  mutate(ss2 = sum(abs(ss)))%>%
  arrange(desc(ss2))%>% 
  pull(gene) %>% unique()

pdf("plots/plot_top_genes_from_PCA.pdf", 
    width= 10, 
    height= 5)
p.list
dev.off()

genes.pc2= enframe(loadings.pca[,2])%>%
  mutate(name= rownames(loadings.pca))%>% 
  arrange(desc(value))

h.man= genes.pc2 %>% 
  top_n(10) %>% pull(name)

top.mm= genes.pc2 %>%
  top_n(10, wt = -value) %>% pull(name)


p1= rank.matrix[h.man,] %>% as.data.frame()%>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "contrast_id", values_to= "logFC")%>%
  left_join(c.df %>% distinct(cc, contrast_id))%>%
  ggplot(., aes(x= cc, y= logFC))+
  facet_grid(~gene)+
  geom_boxplot()+
  geom_jitter(aes(color= contrast_id))+
  geom_hline(yintercept = 0)+ theme_cowplot()

p2= rank.matrix[top.mm,] %>% as.data.frame()%>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "contrast_id", values_to= "logFC")%>%
  left_join(c.df %>% distinct(cc, contrast_id))%>%
  ggplot(., aes(x= cc, y= logFC))+
  facet_grid(~gene)+
  geom_boxplot()+
  geom_jitter(aes(color= contrast_id))+
  geom_hline(yintercept = 0)+theme_cowplot()

cowplot::plot_grid(p1,p2, ncol = 1)

