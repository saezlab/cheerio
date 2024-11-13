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
library(lsa)
library(ComplexHeatmap)
library(circlize)

library(ggrepel)

c.df= readRDS("data/contrasts_query_df_translated3.rds")

c.df$contrast_id%>% unique()%>%
  write.csv("data/contrast_ids.csv")
# number of unique genes

length(unique(c.df$gene))
unique(c.df$contrast_id)

contrast_comparisons<- c.df%>% 
  group_by(gene)%>%
  summarize(contrast_count = n_distinct(contrast_id))

sum(contrast_comparisons$contrast_count >9)

p.hist.counts <- contrast_comparisons%>% 
  ggplot(aes(x= contrast_count))+
  geom_histogram(bins = 30)+
  scale_y_log10()+
  labs(x= "number of contrasts", 
       y= "number of genes")+
  theme_cowplot()

p.hist.counts


# run correlation and hierarchical cluster --------------------------------

get_boxplot_from_matrix<- function(mtx){
  mtx %>% as.data.frame()%>% 
    rownames_to_column("gene")%>%
    pivot_longer(-gene)%>%
    ggplot(aes(x= name, y= value))+
    geom_boxplot()+
    theme_minimal()+
    theme(axis.text.x = element_text(angle= 90, vjust= 0.5, hjust= 1))+
    labs(x= "")
}

unsupervised_wrapper<- function(c.df, 
                                contrast_oi){
  wide_lfc= c.df %>%
    dplyr::select(logFC, contrast_id, gene)%>% 
    filter(contrast_id %in% contrast_oi)%>%
    pivot_wider(names_from = contrast_id, values_from  = logFC, values_fn= mean)%>%
    drop_na() %>%
    as.data.frame()%>% 
    column_to_rownames("gene")%>% 
    as.matrix()
  
  wide_lfc_scaled= (scale(wide_lfc, center = T, scale= T))

  p1<- get_boxplot_from_matrix(wide_lfc)
  
  p2<- get_boxplot_from_matrix(wide_lfc_scaled)
  
  cowplot::plot_grid(p1+labs(y= "log2FC"), p2+labs(y = "scaled log2FC"), align ="h")
  
  ## hierarchical clustering
  
  d= cosine((wide_lfc))
  cl= hclust(as.dist(1-d))
  plot(cl, main = "", xlab = "")
  
  ## correlation
  
  d= cor(wide_lfc, method = "pearson")
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  corr.hmap= ComplexHeatmap::Heatmap(d,row_dend_side = "right", show_column_dend = T,show_row_dend= F, row_names_side = "left",
                                     name = "Pearson's\ncorrelation",
                                     col= col_fun, 
                                     column_names_max_height= unit(20, "cm"),
                                     row_names_max_width= unit(10, "cm"),
                                     #rect_gp = gpar(col = "#303030", lwd = 2,size=2)  # Set the border color
                                     rect_gp = gpar(col = "white", lwd = 1)  # Set the border color
  )
  print(corr.hmap)
  print(dim(wide_lfc))
  
  return(wide_lfc)
}

## define contrasts of interest

# compare all proteomics data:
contrast_prot <- unique(c.df$contrast_id[grepl("prot", c.df$contrast_id)])

#compare all animal models
contrast_mm <- unique(c.df$contrast_id[grepl("mm|rn", c.df$contrast_id)])
mm_res <- unsupervised_wrapper(c.df,contrast_mm )

#get numbers:
cor.test(mm_res[,"Mm_tac_ribo_2d"], mm_res[,"Mm_tac_rna_2d"])
cor.test(mm_res[,"Mm_tac_ribo_2wk"], mm_res[,"Mm_tac_rna_2wk"])
cor.test(mm_res[,"Mm_tac_ribo_2wk"], mm_res[,"Mm_tac_prot_2wk"])
cor.test(mm_res[,"Mm_tac_rna_2wk"], mm_res[,"Mm_tac_prot_2wk"])
cor.test(mm_res[,"Mm_tac_rna_2d"], mm_res[,"Mm_tac_prot_2wk"])

cor.test(mm_res[,"Mm_swim_prot_2wk"], mm_res[,"Mm_tac_prot_2wk"])

#compare all human data
contrast_hs <- unique(c.df$contrast_id[grepl("hs", c.df$contrast_id)])

contrasts_oi <- c.df %>%
  filter(!grepl("HCMvsDCM|snR", contrast_id)) %>%
  filter(grepl("hs", contrast_id)) %>%
  pull(contrast_id)%>% unique()
hs_res <- unsupervised_wrapper(c.df,contrasts_oi )

## cross species
contrast_cs <- c.df %>%
  filter(!grepl("prot|snR", contrast_id)) %>%
  #filter(grepl("mm_TAC", contrast_id)) %>%
  pull(contrast_id)%>% unique()
contrast_cs

cs_res <- unsupervised_wrapper(c.df,contrast_cs )

contrast_cs <- c.df %>%
  filter(!grepl("prot|snR|ribo", contrast_id)) %>%
  #filter(grepl("mm_TAC", contrast_id)) %>%
  pull(contrast_id)%>% unique()
contrast_cs

contrast_cm = c(contrast_cs, "hs_HCMvsNF_snRNA_CM")
cm_res <- unsupervised_wrapper(c.df,contrast_cm )

# PCA  --------------------------------------------------------------------
cm_res_scaled = t(scale((cm_res), center=T, scale = T))

p1<- get_boxplot_from_matrix(cm_res)
p2<- get_boxplot_from_matrix(t(cm_res_scaled))
cowplot::plot_grid(p1+labs(y= "log2FC"), p2+labs(y = "scaled log2FC"), align ="h")

PCA= prcomp((cm_res_scaled), center = F, scale. =F)
#PCA= prcomp(t(cm_res), center = T, scale. =T)

p.df= PCA$x %>% 
  as.data.frame()%>%
  rownames_to_column("c.id")%>%
  as_tibble()%>% 
  mutate(species= ifelse(grepl("hs", c.id), "human", "animal"))
p.df

plot_pca= function(p.df, 
                   pc_x= 1, 
                   pc_y= 2){
  require(ggrepel)
  require(RColorBrewer)
  #RColorBrewer::display.brewer.pal(5, "qual")
  x_col= paste0("PC", pc_x)
  y_col= paste0("PC", pc_y)
  p.df %>% 
  ggplot(aes(x = !!rlang::ensym(x_col),
             y = !!rlang::ensym(y_col), 
             shape= species, 
             color= c.id,
             label = c.id))+
    #geom_text(alpha= 0.6, color="black")+
    geom_text_repel(alpha= 0.6)+
    geom_point(size= 3)+
    theme_cowplot()+
    #scale_color_brewer(type="qual", palette=2)+
    theme(axis.line = element_blank())+
    labs(x= paste0(x_col, " (",as.character(round(PCA$sdev[pc_x]^2/sum(PCA$sdev^2)*100)),"%)"),
         y= paste(y_col, " (",as.character(round(PCA$sdev[pc_y]^2/sum(PCA$sdev^2)*100)),"%)"),
         color= "contrast ID")+
    ggtitle(paste0(""))+
    theme( panel.border = element_rect(colour = "black", fill=NA, size=1))
}

p1= plot_pca(p.df, 1,2)
#p1= plot_pca(p.df, 5,6)
p1

legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12))
  )

p2= plot_pca(p.df, 3,4)
p2

plot_grid(p1+theme(legend.position = "none"),
          p2+theme(legend.position = "none"),
          nrow= 1,
          legend)

p3= plot_pca(p.df, 5,6)

pdf("plots/plot_PCA_contrasts12.pdf",
    width= 8, height= 5)
p1
dev.off()

##screeplot
map(PCA$sdev, function(x){
  round(x^2/sum(PCA$sdev^2)*100, 3)
})%>% unlist()%>% plot()


## for CM comparison. pathologic hypertrophy on PC1 (negative) 
## and species on PC2 negative is mouse
## Q which are the genes that are species independent on pathologic hypertrophy
hist(PCA$rotation[,1], breaks = 50)
hist(PCA$rotation[,2], breaks = 50)


p.pcaloadings= PCA$rotation[,1:5] %>% 
  as_data_frame() %>% 
  mutate(gene= rownames(PCA$rotation))

p.pcaloadings%>%
  ggplot(., aes(x= PC3, y= PC4))+
  geom_point()
  
p.pcaloadings<-
p.pcaloadings %>% 
  mutate(mm= (PC1 < -0.03 ) & (PC2 < 0.03),
         hs= (PC1 < -0.03 ) & (PC2 > 0.03), 
         label= ifelse(hs | mm, gene, ""),
         prio= -(PC1)*abs(PC2),
         prio2= -PC1* (1-PC2))
plot.loads= p.pcaloadings%>%
  ggplot(., aes(x= PC1, y= PC2, color= prio))+
  geom_label_repel(aes(label= label),
                   alpha= 0.6,
                   color= "black",
                   max.overlaps = 40)+
  geom_point()+#color= "darkgrey")+
  theme_cowplot()+
  scale_color_gradient2(low= "black",mid= "darkgrey", high = "red")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=1))+
  coord_equal()+
  labs(color= "")
plot.loads

pdf("plots/plot_PCA_contrasts12_loadings.pdf",
    width= 8, height= 5)
plot.loads
dev.off()

#get genes that associate with hypertrophy but not w species

p.pcaloadings%>%
  arrange(desc(prio))%>%
  slice_head(n= 20)%>% pull(gene)

# pc1.df= enframe(sort(PCA$rotation[,1:2]))%>% arrange(desc(abs(value)))
# 
# species_diff %>% 
#   mutate(dir= sign(value))%>%
#   group_by(dir)%>%
#   mutate(value2= abs(value)) %>% 
#   top_n(n = 20,wt = value2)%>%
#   arrange(desc(value))%>% 
#   print(n=100)



# NMF-------------------------------------------------------------------------
library(NMF)
library(ComplexHeatmap)
data <- as.matrix(r.m) + 11
get_boxplot_from_matrix(data)
estimate <- nmfEstimateRank(data, range = 3:6)
set.seed(1)
plot(estimate)
estimate$consensus$`4`
nmf_result <- nmf(data, rank = 6, method = 'brunet', nrun = 10)
W <- basis(nmf_result)
H <- coef(nmf_result)

Heatmap(t(H), cluster_columns = F)

sort(W[,2], decreasing = T)[1:30]

get_boxplot_from_matrix(t(r.m["ANKRD1",]))
plot(W["ANKRD1", ])


library(msigdbr)
df<- msigdbr::msigdbr(category = "H")%>%
  distinct(gs_name, gene_symbol)

res<- decoupleR::run_ulm(mat= W, network = df, .source="gs_name", .target = "gene_symbol")
pways <- res %>% filter(p_value< 0.0001)%>% pull(source)

res %>% 
  filter(source %in% pways) %>% 
  ggplot(aes(x= condition, y= source, fill =score))+
  geom_tile()+
  scale_fill_gradient2()

# Visualize the correlation matrix
heatmap(cor_matrix)
# plot correlation of logFC -----------------------------------------------

corrplot::corrplot(corr = cor(r.m, method = "pearson"))

# add hierachical cluster with cosine

d= cosine((rank.matrix))
cl= hclust(as.dist(1-d))
plot(cl, main = "", xlab = "")

pdf("plots/plot_clust_dend_w_prot.pdf",
    height= 6,
    width= 4)
plot(cl, main = "", xlab = "")
dev.off()

d= cor(rank.matrix, method = "spearman")
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

corr.hmap= ComplexHeatmap::Heatmap(d,row_dend_side = "right", show_column_dend = T,show_row_dend= F, row_names_side = "left",
                        name = "spearman rho",
                        col= col_fun, 
                        column_names_max_height= unit(20, "cm"),
                        row_names_max_width= unit(10, "cm"),
                        #rect_gp = gpar(col = "#303030", lwd = 2,size=2)  # Set the border color
                        rect_gp = gpar(col = "white", lwd = 1)  # Set the border color
                      )
corr.hmap
pdf("plots/plot_corr_matrix.pdf", 
    width= 6.5, height= 6)
corr.hmap
dev.off()



# wgcna -------------------------------------------------------------------
library(WGCNA)

logFC_matrix<- r.m

powers = c(1:20)

# Pick the best soft-threshold power based on scale-free topology fit
sft = pickSoftThreshold(logFC_matrix, powerVector = powers, verbose = 5)

# Plot the results to visually inspect the optimal power
par(mfrow = c(1, 1))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     type = "n", xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
softPower = sft$powerEstimate  # Chosen power based on previous analysis
softPower = 14
# Construct the adjacency matrix
adjacency = adjacency(logFC_matrix, power = softPower)

# Turn the adjacency matrix into a TOM
TOM = TOMsimilarity(adjacency)

# Turn the TOM into a dissimilarity matrix
dissTOM = 1 - TOM

# Perform hierarchical clustering
geneTree = hclust(as.dist(dissTOM), method = "average")
plot(geneTree)
# Dynamic tree cut to identify modules
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            cutHeight = 1.3, 
                            #deepSplit = 2,
                            #pamRespectsDendro = FALSE,
                            minClusterSize = 10)

# Convert numeric labels into colors for visualization
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

#Plot the dendrogram and module colors
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
# Calculate the module eigengenes (first principal component of each module)
MEs = moduleEigengenes(logFC_matrix, colors = dynamicColors)$eigengenes

# Correlate each module eigengene with contrasts
moduleTraitCor = cor(MEs, logFC_matrix, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(logFC_matrix))

# Visualize the correlations as a heatmap
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(logFC_matrix),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = ""),
               setStdMargins = FALSE,
               cex.text = 0.7, zlim = c(-1,1))
# interpret pc2 -----------------------------------------------------------

msigDB= readRDS("/home/jan/R-projects/scell_hfpef/data/prior_knowledge/Genesets_Dec19.rds")
library(decoupleR)
PCA$center
loadings.pca= PCA$rotation

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

msig_res= decoupleR::run_ulm(mat = loadings.pca[,1:2], net= msigDBnet)

write.csv(msig_res, "data/misg_enriched_in_pca_loadings.csv")

pways= msig_res %>%filter(statistic=="ulm") %>% mutate(p_adj= p.adjust(p_value))%>%
  group_by(condition, source)%>%
  filter(any(p_adj<0.001))%>%
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



# plot n# genes -----------------------------------------------------------

p.coverag <- c.df%>% 
  distinct(contrast_id, max.ranks)%>%
  ggplot(., aes(x= reorder(contrast_id,max.ranks), y= max.ranks))+
  geom_col()+
  coord_flip()+
  theme_cowplot()+
  labs(y= "gene coverage", 
       x= "")

p.deg <- c.df %>%
  count(sig) %>%
  filter(sig)%>%
  ggplot(., aes(x= reorder(contrast_id,n), y= n))+
  geom_col()+
  coord_flip()+
  theme_minimal()+
  labs(y= "number of DEGs", 
       x= "")

pdf("plots/plot_coverage.pdf",
    width = 15, 
    height= 10)
plot_grid(p.coverag, p.deg)
dev.off()
  


# get signature -----------------------------------------------------------

# generate signatures -----------------------------------------------------


library(survcomp)
# update the 
c.df2<- c.df%>% 
  group_by(contrast_id)%>%
  mutate(pval_min = min(pval[pval > 0]),  ## modify FDR values that are INF to smallest non-zero FDR within that contrast
         pval_mod= ifelse(pval== 0, pval_min, pval))

contrasts_oi<- contrast_cm
missing_na<- 2
(length(contrasts_oi)-missing_na)
x<- c.df2 %>% 
  filter(contrast_id %in% contrasts_oi)%>%
  group_by(gene)%>% 
  select(gene, pval_mod)%>%
  mutate(gene_count = as.numeric(ave(gene, gene, FUN = length)))%>%
  filter(gene_count > (length(contrasts_oi)-missing_na) )%>%
  mutate(fisher_p=  survcomp::combine.test(pval_mod, "fisher", na.rm = T))%>%
  distinct(gene, fisher_p)%>%
  ungroup()%>% 
  mutate(fisher_p_adj= p.adjust(fisher_p))


ggplot(x, aes(x= reorder(gene, fisher_p), y= -log10(fisher_p)))+
  geom_point()+
  theme_minimal()+
  theme(axis.text.x=element_blank())



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

x= get_top_consistent_gene(joint_contrast_df = c.df, 
                           query_contrasts =contrasts_oi, 
                           missing_prop = 20 )
x

