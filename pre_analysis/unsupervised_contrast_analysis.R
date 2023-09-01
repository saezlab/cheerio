## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-06-14
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## breakdown programs


library(tidyverse)
library(cowplot)
library(ggrepel)

c.df= readRDS("data/contrasts_query_df.rds")

unique(c.df$contrast_id)
contrast_focus= c("Mm_tac_rna_2wk", 
                  "Mm_swim_rna_2wk", 
                  "Mm_swim_rna_2d", 
                  "Mm_tac_rna_2d", "Hs_fetal_Akat14", "Hs_ReHeaT",
                  "Hs_fetal_Spurell22", 
                  #"Rn_invitro_ribo",
                  "Rn_invitro_rna", 
                  "Hs_singlecell_HCMvsNF_Cardiomyocyte",
                  #"Hs_singlecell_HCMvsNF_Endothelial I",
                 # "Hs_singlecell_HCMvsNF_Macrophage",
                  "Hs_bulk_HCMvsNF"
                  #"Hs_singlecell_HCMvsNF_Fibroblast"
                 )

contrast_focus= c("Mm_tac_rna_2wk", 
                  "Mm_swim_rna_2wk", 
                  "Mm_swim_rna_2d", 
                  "Mm_tac_rna_2d", "Hs_fetal_Akat14", "Hs_ReHeaT",
                  "Hs_fetal_Spurell22", 
                  #"Rn_invitro_ribo",
                  #"Rn_invitro_rna", 
                  #"Hs_singlecell_HCMvsNF_Cardiomyocyte",
                  #"Hs_singlecell_HCMvsNF_Endothelial I",
                  # "Hs_singlecell_HCMvsNF_Macrophage",
                  "Hs_bulk_HCMvsNF",
                  "Hs_singlecell_HCMvsNF_Fibroblast"
)


contrast_focus= c("Mm_tac_rna_2wk", 
                  "Mm_swim_rna_2wk", 
                  "Mm_tac_rna_2d",
                  "Mm_swim_rna_2d",
                  "Rn_invitro_ribo")#"Rn_invitro_rna")

rank.df= c.df %>%
  mutate(min.FDR = min(FDR[FDR > 0]),
         FDR_mod= ifelse(FDR== 0, min.FDR, FDR), 
         new_weight= (logFC)*-log10(FDR_mod))%>%
  #filter(FDR<.05)%>%
  select(logFC, contrast_id, gene)%>% 
  filter(contrast_id %in% contrast_focus)%>%
  pivot_wider(names_from = contrast_id, values_from  = logFC, values_fn= mean)

rank.df


                       
rank.matrix=rank.df %>% 
  drop_na() %>%
  #filter(!is.na(gene))%>%
  as.data.frame()%>% 
  column_to_rownames("gene")%>% as.matrix()

#ceck distribution of gex values
#apply(rank.matrix,2, hist)

#rank.matrix
dim(rank.matrix)
r.m= scale(rank.matrix)

#apply(r.m, 2, hist)

PCA= prcomp(t(r.m), center = F, scale. =F)

p.df= PCA$x %>% 
  as.data.frame()%>%
  rownames_to_column("c.id")%>%
  as_tibble()

plot_pca= function(p.df, 
                   pc_x, 
                   pc_y){
  
  x_col= paste0("PC", pc_x)
  y_col= paste0("PC", pc_y)
  p.df %>% 
  ggplot(aes(x = !!rlang::ensym(x_col), y = !!rlang::ensym(y_col), color= c.id))+
    geom_point(size= 3) +
    theme_minimal()+
    labs(x= paste0(x_col, " (",as.character(round(PCA$sdev[pc_x]^2/sum(PCA$sdev^2)*100)),"%)"),
         y= paste(y_col, " (",as.character(round(PCA$sdev[pc_y]^2/sum(PCA$sdev^2)*100)),"%)"))+
    ggtitle(paste0(""))+
    geom_label_repel(aes(label= c.id))+
    theme_cowplot()
}

p1= plot_pca(p.df, 1,2)
p2= plot_pca(p.df, 3,4)
p3= plot_pca(p.df, 5,6)

plot_grid(p1,p2,p3)

map(PCA$sdev, function(x){
  round(x^2/sum(PCA$sdev^2)*100, 3)
})%>% unlist()%>% plot()


#PC1 associates with species

hist(PCA$rotation[,1], breaks = 50)

plot(rank.matrix["MYH10",])
rank.matrix["STAT3",]


## for CM comparison. pathologic hypertrophy on PC1 (negative) 
## and species on PC2 negative is mouse
## Q which are the genes that are species independent on pathologic hypertrophy
hist(PCA$rotation[,1], breaks = 20)
hist(PCA$rotation[,2], breaks = 20)


p.pcaloadings= PCA$rotation[,1:5] %>% 
  as_data_frame() %>% 
  mutate(gene= rownames(PCA$rotation))

p.pcaloadings%>%
  ggplot(., aes(x= PC3, y= PC4))+
  geom_point()
  
  
p.pcaloadings %>% 
  mutate(mm= (PC1 < -0.03 ) & (PC2 < 0.03),
         hs= (PC1 < -0.03 ) & (PC2 > 0.03), 
         label= ifelse(hs | mm, gene, ""),
         prio= -(PC1)*abs(PC2))%>%
  #pull(prio)%>% hist()

  ggplot(., aes(x= PC1, y= PC2, color= prio))+
  geom_point()+
  geom_label_repel(aes(label= label))+
  theme_cowplot()+
  scale_color_gradient2(low= "black",mid= "darkgrey", high = "red")

r.m["XIRP2",] %>% enframe(name = "contrast_id")%>%
  left_join(c.df %>% distinct(cc, contrast_id))%>%
  ggplot(., aes(x= cc, y= value))+
  geom_boxplot()+
  geom_jitter(aes(color= contrast_id))
r.m["PHLDA1",]



# pc1.df= enframe(sort(PCA$rotation[,1:2]))%>% arrange(desc(abs(value)))
# 
# species_diff %>% 
#   mutate(dir= sign(value))%>%
#   group_by(dir)%>%
#   mutate(value2= abs(value)) %>% 
#   top_n(n = 20,wt = value2)%>%
#   arrange(desc(value))%>% 
#   print(n=100)


# interpret pc2 -----------------------------------------------------------

msigDB= readRDS("/home/jan/R-projects/sc_hfpef/data/prior_knowledge/Genesets_Dec19.rds")
library(decoupleR)

loadings.pca= p.pcaloadings%>% as.data.frame()%>%
  column_to_rownames("gene")

msig_sel = msigDB[c("MSIGDB_HMARKS",
         "MSIGDB_REACTOME",
         #"MSIGDB_TF" ,
         "MSIGDB_CANONICAL",
         "MSIGDB_KEGG",
         "MSIGDB_BIOCARTA" )]
msigDBnet= lapply(names(msig_sel), function(y){
  head(y)
  enframe(msig_sel[[y]], name="source", value= "target")%>% 
    #mutate(collection= y) %>%
    unnest(target)%>% mutate(mor=1)
  })%>% do.call(rbind, .)%>% distinct(source, target,mor)

msig_res= decoupleR::run_ulm(mat = loadings.pca, net= msigDBnet)

pways= msig_res %>%filter(statistic=="ulm") %>% mutate(p_adj= p.adjust(p_value))%>%
  group_by(condition, source)%>%
  filter(any(p_adj<0.05))%>%
  pull(source
       )

px= msig_res %>% mutate(p_adj= p.adjust(p_value)) %>%
  filter(source %in% pways)%>%
  ggplot(aes(x=score, y= reorder(source,score), fill = condition))+
  facet_grid(~ condition)+
  geom_col(aes(color= p_adj<0.01))+
  scale_color_manual(values= c("TRUE"= "black", "FALSE"= "white"))+
  theme_cowplot()+
  geom_vline(xintercept = 0)

pdf("big_pc_hmap.pdf",width= 20, height= 25)
px
dev.off()


msig_res %>% mutate(p_adj= p.adjust(p_value)) %>%
  filter(source %in% pways, 
         condition== "PC2",
         score>0)%>% print(n=100)

genes.pc2= enframe(loadings.pca[,2])%>%mutate(name= rownames(loadings.pca))%>% 
  arrange(desc(value))

h.man= genes.pc2 %>% 
   top_n(10) %>% pull(name)

top.mm= genes.pc2 %>%
  top_n(10, wt = -value) %>% pull(name)



p1= r.m[h.man,] %>% as.data.frame()%>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "contrast_id", values_to= "logFC")%>%
  left_join(c.df %>% distinct(cc, contrast_id))%>%
  ggplot(., aes(x= cc, y= logFC))+
  facet_grid(~gene)+
  geom_boxplot()+
  geom_jitter(aes(color= contrast_id))+
  geom_hline(yintercept = 0)+ theme_cowplot()

p2= r.m[top.mm,] %>% as.data.frame()%>% rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "contrast_id", values_to= "logFC")%>%
  left_join(c.df %>% distinct(cc, contrast_id))%>%
  ggplot(., aes(x= cc, y= logFC))+
  facet_grid(~gene)+
  geom_boxplot()+
  geom_jitter(aes(color= contrast_id))+
  geom_hline(yintercept = 0)+theme_cowplot()

cowplot::plot_grid(p1,p2, ncol = 1)
# plot n# genes -----------------------------------------------------------

c.df%>% 
  distinct(contrast_id, max.ranks)%>%
  ggplot(., aes(x= reorder(contrast_id,max.ranks), y= max.ranks))+
  geom_col()+
  coord_flip()+
  theme_cowplot()

c.df %>%
  count(sig) %>%
  filter(sig)%>%
  ggplot(., aes(x= reorder(contrast_id,n), y= n))+
  geom_col()+
  coord_flip()+
  theme_cowplot()

library(cowplot)



# plot correlation of logFC -----------------------------------------------

rank.df = c.df %>%
  mutate(min.FDR = min(FDR[FDR > 0]),
         FDR_mod= ifelse(FDR== 0, min.FDR, FDR), 
         new_weight= sign(logFC)*-log10(FDR_mod))%>%
  #filter(FDR<.05)%>%
  select(logFC, contrast_id, gene)%>% 
  filter(contrast_id %in% contrast_focus)%>%
  pivot_wider(names_from = contrast_id, values_from  = logFC, values_fn= mean)

rank.matrix=rank.df %>% 
  drop_na() %>%
  #filter(!is.na(gene))%>%
  as.data.frame()%>% 
  column_to_rownames("gene")%>% as.matrix()


corrplot::corrplot(corr = cor(r.m, method = "pearson"))

# add hierachical cluster with cosine
library(lsa)
sc.m= scale((rank.matrix))
d= cosine(r.m)
cl= hclust( as.dist(1-d))
plot(cl)

d1= dist(t(sc.m))
plot(hclust(d1))

ComplexHeatmap::Heatmap(d)

dim(r.m)
