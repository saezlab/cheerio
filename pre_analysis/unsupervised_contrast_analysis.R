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
                  "Rn_invitro_ribo",
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
  select(new_weight, contrast_id, gene)%>% 
  filter(contrast_id %in% contrast_focus)%>%
  pivot_wider(names_from = contrast_id, values_from  = new_weight, values_fn= mean)

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
r.m= apply(rank.matrix, 2, scale)
rownames(r.m)= rownames(rank.matrix)
r.m= scale(rank.matrix, center = T, scale= T)
#apply(r.m, 2, hist)

PCA= prcomp(t(r.m), center = F, scale. =F)

p.df= PCA$x %>% 
  as.data.frame()%>%
  rownames_to_column("c.id")%>%
  as_tibble()


p.df %>% 
  ggplot(., aes(x= PC1, y= PC2, color= c.id))+
  geom_point(size= 3) +
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))+
  geom_label_repel(aes(label= c.id))+
  theme_cowplot()

p.df %>% 
  ggplot(., aes(x= PC3, y= PC4, color= c.id))+
  geom_point(size= 3) +
  theme_minimal()+
  labs(x= paste0("PC3 (",as.character(round(PCA$sdev[3]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC4 (",as.character(round(PCA$sdev[4]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))+
  geom_label_repel(aes(label= c.id))+
  theme_cowplot()

map(PCA$sdev, function(x){
  round(x^2/sum(PCA$sdev^2)*100, 3)
})%>% unlist()%>% plot()


#PC1 associates with species

hist(PCA$rotation[,1], breaks = 20)

plot(rank.matrix["MYH10",])
rank.matrix["STAT3",]


## for CM comparison. pathologic hypertrophy on PC1 (negative) 
## and species on PC2 negative is mouse
## Q which are the genes that are species independent on pathologic hypertrophy
hist(PCA$rotation[,1], breaks = 20)
hist(PCA$rotation[,2], breaks = 20)

PCA$rotation
p.pcaloadingS= PCA$rotation[,1:4] %>% 
  as_data_frame() %>% 
  mutate(gene= rownames(PCA$rotation))

p.pcaloadingS%>%
  ggplot(., aes(x= PC3, y= PC4))+
  geom_point()
  geom_te

  
p.pcaloadingS %>% 
  mutate(mm= (PC1< -0.05 ) & (PC2< -0.1),
         hs= (PC1< -0.05 ) & (PC2> -0.05), 
         label= ifelse(hs | mm, gene, ""))%>%
  ggplot(., aes(x= PC1, y= PC2))+
  geom_point()+
  geom_label_repel(aes(label= label))

r.m["NPPB",] %>% enframe(name = "contrast_id")%>%
  left_join(c.df %>% distinct(cc, contrast_id))%>%
  ggplot(., aes(x= cc, y= value))+
  geom_boxplot()+
  geom_jitter(aes(color= contrast_id))
r.m["PHLDA1",]

hist(r.m[,8])
filter(mm | hs)
pc1.df= enframe(sort(PCA$rotation[,1:2]))%>% arrange(desc(abs(value)))

species_diff %>% 
  mutate(dir= sign(value))%>%
  group_by(dir)%>%
  mutate(value2= abs(value)) %>% 
  top_n(n = 20,wt = value2)%>%
  arrange(desc(value))%>% 
  print(n=100)

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
