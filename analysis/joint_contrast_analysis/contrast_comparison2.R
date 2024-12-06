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
library(ggdendro)

c.df= readRDS("app_data/contrasts_query_df_translated3.rds")
source("sub/global.R")

c.df$contrast_id%>% unique()%>%
  write.csv("data/contrast_ids.csv")
# number of unique genes

length(unique(c.df$gene))
unique(c.df$contrast_id)

cc_colors <- c("A" = "#20558c",
               "B" = "#ff8708",
               "C" = "#a8629d",
               "D" = "#0c9659")
#aesthetics and meta data
label_colors <- c.df %>%
  distinct(contrast_id, cc, model, modal) %>%
  arrange(cc, contrast_id)%>%
  mutate(color = cc_colors[as.character(cc)], 
         species= str_to_title(substr(contrast_id, 1, 2)),
         hypertrophy= "pathologic",
         hypertrophy= ifelse(cc %in% c("C", "D"), NA, hypertrophy),
         hypertrophy = ifelse(grepl("swim",contrast_id), "physiologic", hypertrophy),
  )


# run correlation and hierarchical cluster --------------------------------

get_boxplot_from_matrix<- function(mtx){
  mtx %>% as.data.frame()%>% 
    rownames_to_column("gene")%>%
    pivot_longer(-gene)%>%
    ggplot(aes(x= name, y= value))+
    geom_boxplot()+
    theme_cowplot()+
    theme(axis.text.x = element_text(angle= 90, vjust= 0.5, hjust= 1))+
    labs(x= "")
}
contrast_oi<- contrast_cs

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
  #plot(cl, main = "", xlab = "")
  
  # Perform hierarchical clustering
  hc <- hclust(as.dist(1 - d))
  
  # Convert clustering to a dendrogram object
  dendro_data <- ggdendro::dendro_data(as.dendrogram(hc))
  
  # Plot using ggplot2
  p.dend<- ggplot(segment(dendro_data)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_text(data = label(dendro_data),
              aes(x = x, y = -0.05, label = label), # y = 0 to align labels at the same baseline
              hjust = 1, angle = 90, size = 3.5) + # Customize text size and angle
    theme_bw() +
    labs(title = "",
         x = "",
         y = "Height") +
    theme(axis.text.x = element_blank(), # Remove x-axis text
          axis.ticks = element_blank(),
          plot.margin = margin(20, 20, 120, 20, "pt"),
          panel.border = element_blank())+
    coord_cartesian(clip = 'off') # Remove x-axis ticks
  
  print(p.dend)
  ## correlation
  #calc the corr
  d= cor(wide_lfc, method = "pearson")
  hc2 <- hclust(as.dist(1 - d))
  plot(hc2)
  # Custom distance function for clustering
  dist_function = function(x) as.dist(1 - x)
  
  # aesthetics
  # filling
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  # labels
  
  modal_colors <- c("transcriptome" = "#FFA07A", "proteome" = "#20B2AA", "translatome" = "#FFD700") # 3 levels
  model_colors <- c("pathologic" = "#20558c", "physiologic" = "#FFD700", "NA"= "grey")   # 2 levels + NA
  species_colors <- c("Hs" = "#ff8708", "Rn" = "#d62f2f", "Mm" = "#0c9659")        # 3 species
  
  
  # Create named vectors for each annotation variable
  contrast_colors <- setNames(label_colors$color, label_colors$contrast_id)
  
  row_colors <- contrast_colors[rownames(d)]
  column_colors <- contrast_colors[colnames(d)]
  
  column_cc_vector <- setNames(label_colors$cc, label_colors$contrast_id)
  column_modal_vector <- setNames(label_colors$modal, label_colors$contrast_id)
  column_model_vector <- setNames(label_colors$hypertrophy, label_colors$contrast_id)
  column_species_vector <- setNames(label_colors$species, label_colors$contrast_id)
  
  # Match the order of the column names in the heatmap matrix
  column_cc_data <- column_cc_vector[colnames(d)]
  column_modal_data <- column_modal_vector[colnames(d)]
  column_model_data <- column_model_vector[colnames(d)]
  column_species_data <- column_species_vector[colnames(d)]
  
  # Define the column annotation
  column_annotation <- HeatmapAnnotation(
    Category = column_cc_data,       # CC annotation
    Modality = column_modal_data, # Modal annotation
    Hypertrophy = column_model_data, # Model annotation
    Species = column_species_data, # Species annotation
    col = list(
      Category = cc_colors,
      Modality = modal_colors,
      Hypertrophy = model_colors,
      Species = species_colors
    ),
    show_legend = F
    # annotation_legend_param = list(
    #   title_gp = gpar(fontsize = 10, fontface = "bold"),  # Styling for legend titles
    #   grid_height = unit(5, "mm")  # Height of the legend grids
    # )
  )
  corr.hmap= ComplexHeatmap::Heatmap(d,
                                     row_dend_side = "right",
                                     show_column_dend = T,
                                     show_row_dend= F, 
                                     row_names_side = "left",
                                     name = "Pearson's\ncorrelation",
                                     col= col_fun, 
                                     cluster_rows = hclust(dist_function(d)),  # Custom row clustering
                                     cluster_columns = hclust(dist_function(t(d))),
                                     clustering_method_columns = "average", 
                                     clustering_method_rows =  "average",
                                     column_names_max_height= unit(20, "cm"),
                                     row_names_max_width= unit(10, "cm"),
                                     column_dend_height = unit(40, "mm"),
                                     column_names_side = "top",
                                     border = TRUE,
                                     show_heatmap_legend = F,
                                     border_gp = gpar(col = "black"),
                                     #rect_gp = gpar(col = "#303030", lwd = 2,size=2)  # Set the border color
                                     rect_gp = gpar(col = "white", lwd = 1),   # Set the border color
                                     #row_names_gp = gpar(col = row_colors),    # Apply colors to row names
                                     #column_names_gp = gpar(col = column_colors), # Apply colors to column names
                                     ##annotations:
                                     top_annotation = column_annotation 
                                     
  )
  
  print(corr.hmap)
  print(dim(wide_lfc))
  # Pearson correlation legend
  heatmap_legend = Legend(
    title = "Pearson's\ncorrelation",
    col_fun = col_fun,
    grid_width = unit(5, "mm"),
    grid_height = unit(20, "mm")
  )
  
  # CC category legend
  cc_legend = Legend(
    labels = names(cc_colors),  # The labels for the legend
    legend_gp = gpar(fill = cc_colors),  # Colors for the legend dots
    title = "Data\ncategory",
    grid_width = unit(5, "mm"),
    grid_height = unit(5, "mm")
  )
  
  # Modal legend
  modal_legend = Legend(
    labels = names(modal_colors),
    legend_gp = gpar(fill = modal_colors),
    title = "Modality",
    grid_width = unit(5, "mm"),
    grid_height = unit(5, "mm")
  )
  
  # Model legend
  model_legend = Legend(
    labels = names(model_colors),
    legend_gp = gpar(fill = model_colors),
    title = "Hypertrophy",
    grid_width = unit(5, "mm"),
    grid_height = unit(5, "mm")
  )
  
  # Species legend
  species_legend = Legend(
    labels = names(species_colors),
    legend_gp = gpar(fill = species_colors),
    title = "Species",
    grid_width = unit(5, "mm"),
    grid_height = unit(5, "mm")
  )
  
  # Combine all legends
  combined_legend = packLegend(
    heatmap_legend,
    cc_legend,
    modal_legend,
    model_legend,
    species_legend,
    direction = "vertical",  # Arrange legends vertically
    gap = unit(5, "mm")      # Space between legends
  )
  # Draw the heatmap with adjusted legend placement
  p.hmap<- draw(corr.hmap, annotation_legend_list = combined_legend, 
                annotation_legend_side = "left" )
  # Draw the heatmap with the stacked legends
  
  
  return(list("lfc_mtx"= wide_lfc, 
              "p.den"= p.dend, 
              "p.hmap"= p.hmap,
              "corr.matrix"= d)
  )
}

## define contrasts of interest

## MAIN CONTRAST

#compare all animal models
contrast_mm <- unique(c.df$contrast_id[grepl("mm|rn", c.df$contrast_id)])

## cross species
contrast_cs <- c.df %>%
  #filter(!grepl("prot|snR", contrast_id)) %>%
  filter(!grepl("snR", contrast_id) | grepl("HCMvsNF_snRNA_CM", contrast_id)) %>%
  #filter(!grepl("HCMvsDCM", contrast_id)) %>%
  pull(contrast_id)%>% unique()
contrast_cs

cs_res <- unsupervised_wrapper(c.df,contrast_cs )

#save correlation results:
cs_res$corr.matrix%>% 
  as.data.frame()%>%
  rownames_to_column("comparison")%>%
  write_csv("analysis/joint_contrast_analysis/correlation_matrix_all.csv")
dim(cs_res$lfc_mtx)

pdf("figures/corr_hmap_cross_species.pdf",
    width= 8.2, height= 8.5)
cs_res$p.hmap
dev.off()

## no pro
contrast_cm <- c.df %>%
  #filter(!grepl("prot|snR|ribo|HFvsNF|HCMvsDCM", contrast_id)) %>%
  filter(!grepl("prot|snR|HCMvsDCM", contrast_id)) %>%
  #filter(grepl("mm_TAC", contrast_id)) %>%
  pull(contrast_id)%>% unique()
contrast_cm
contrast_cm = c(contrast_cm, "hs_HCMvsNF_snRNA_CM")
cm_res <- unsupervised_wrapper(c.df,contrast_cm )

cm_res$corr.matrix%>% 
  as.data.frame()%>%
  rownames_to_column("comparison")%>%
  write_csv("analysis/joint_contrast_analysis/correlation_matrix_noprot.csv")

pdf("figures/corr_hmap_noprot.pdf",
    width= 8, height= 8)
cm_res$p.hmap
dev.off()
## save dendogram


pdf("figures/dendogram_cosine.pdf",
    width= 8, height= 4)
cowplot::plot_grid(cs_res$p.den, cm_res$p.den,
                   ncol= 2, rel_widths= c(1,0.7))
dev.off()


## only patho
contrast_patho<- contrast_cm[!grepl( "swim", contrast_cm)]

cm_res2 <- unsupervised_wrapper(c.df,contrast_patho )
cm_res2$p.hmap
cm_res2$corr.matrix%>% 
  as.data.frame()%>%
  rownames_to_column("comparison")%>%
  write_csv("analysis/joint_contrast_analysis/correlation_matrix_onlypatho.csv")
# PCA  --------------------------------------------------------------------
plot_pca= function(p.df,
                   PCAres, 
                   pc_x= 1, 
                   pc_y= 2){
  require(ggrepel)
  #require(RColorBrewer)
  
  loadings <- as.data.frame(PCAres$rotation) %>%
    rownames_to_column("gene") %>%
    filter(gene %in% c("NPPA", "NPPB", 
                       "MYH7", "ATP2A2"
                       ))%>%
    mutate(across(where(is.numeric), ~ .x * 1000))
  
  
  #RColorBrewer::display.brewer.pal(5, "qual")
  x_col= paste0("PC", pc_x)
  y_col= paste0("PC", pc_y)
  my_palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  
  p <- p.df %>% 
    mutate(hypertrophy= factor(ifelse(is.na(hypertrophy), 
                               "NA", 
                               hypertrophy),
                               levels= c("physiologic", 
                                         "pathologic",
                                         "NA"))
           )%>%
    ggplot(aes(x = !!rlang::ensym(x_col),
               y = !!rlang::ensym(y_col), 
               shape= hypertrophy, 
               color= species,
               label = contrast_id))+
    #geom_text(alpha= 0.6, color="black")+
    geom_text_repel(alpha= 1, 
                    min.segment.length = unit(1,"mm"),
                     force_pull = 0.01, 
                    force =20,
                     size= 3, 
                     color="black")+
    geom_point(size= 3)+
    theme_cowplot()+
    scale_color_manual(values=  c("Hs" = "#ff8708", 
                                  "Rn" = "#d62f2f", 
                                  "Mm" = "#0c9659"))+
    theme(axis.line = element_blank())+
    labs(x= paste0(x_col, " (",as.character(round(PCAres$sdev[pc_x]^2/sum(PCAres$sdev^2)*100)),"%)"),
         y= paste(y_col, " (",as.character(round(PCAres$sdev[pc_y]^2/sum(PCAres$sdev^2)*100)),"%)"),
         color= "Species",
         shape= "Hypertrophy")+
    ggtitle(paste0(""))+
    theme( panel.border = element_rect(colour = "black", fill=NA, size=1))
  
  p <- p +
    geom_segment(data = loadings, 
                 aes(x = deframe(p.df[p.df$contrast_id=="ZeroLFC", "PC1"]),
                     y = deframe(p.df[p.df$contrast_id=="ZeroLFC", "PC2"]),
                     xend = PC1, yend = PC2), # Scale for visibility
                 arrow = arrow(length = unit(0.2, "cm")), 
                 color = "darkgrey",
                 inherit.aes = FALSE ) +
    geom_text_repel(data = loadings, 
                    aes(x = PC1, y = PC2, label = gene), 
                    color = "darkgrey", inherit.aes = FALSE ,
                    size = 3)
  p
}

pca_wrap<- function(c_output){
  cm_res_scaled = t(scale((c_output$lfc_mtx), center=T, scale = T))
  dim(cm_res_scaled)
  p1<- get_boxplot_from_matrix(c_output$lfc_mtx)
  p2<- get_boxplot_from_matrix(t(cm_res_scaled))
  p3<- get_boxplot_from_matrix(scale(t(cm_res_scaled), center=T, scale=T))
  lfcWithZero <- rbind(cm_res_scaled, ZeroLFC=0)
  dim(lfcWithZero)
 
  p.lfc<-cowplot::plot_grid(p1+labs(y= "log2FC"), p2+labs(y = "scaled log2FC"), align ="h")
  p.lfc
  
  
  PCA= prcomp((lfcWithZero), center = T, scale. =T)
  #PCA= prcomp(t(cm_res), center = T, scale. =T)
  
  p.df= PCA$x %>% 
    as.data.frame()%>%
    rownames_to_column("contrast_id")%>%
    as_tibble()%>% 
    #mutate(species= ifelse(grepl("hs", contrast_id), "human", "animal"))%>%
    left_join(label_colors, by= "contrast_id")
  

  p1= plot_pca(p.df = p.df,PCA, 1,2)
  p2= plot_pca(p.df,PCA, 3,4)
  
  list(PCA, 
       p.lfc,
       p.df, 
       p1+ ggtitle(dim(cm_res_scaled)[2]))
  #p2)
}

pca_res<- pca_wrap(c_output = cm_res)
pca_res

pdf("figures/pca_scale_logfc.pdf",
    width= 6, height= 5)
pca_res[[2]]
dev.off()

pdf("figures/pca_biplot.pdf",
    width= 6, height= 6)
pca_res[[4]]+coord_equal()
dev.off()


##screeplot
expl_var<-map(pca_res[[1]]$sdev, function(x){
  round(x^2/sum(PCA$sdev^2)*100, 3)
})%>% unlist()

names(expl_var)= paste0("PC ", 1:length(expl_var))
p.scree<-enframe(expl_var)%>%
  mutate(name= factor(name, levels=paste0("PC ", 1:length(expl_var))))%>%
  ggplot(aes(x= name, y= value))+
  geom_col(width = 0.7, 
           fill ="grey", 
           color= "black")+
  theme_cowplot()+
  labs(x="", y= "Explained variance (%)")+
  theme(axis.text.x = element_text(angle= 90, hjust = 0, vjust = 0.5))

p.scree
pdf("figures/pca_screeplot.pdf",
    width = 3, height= 3)
p.scree
dev.off()

## for CM comparison. pathologic hypertrophy on PC1 (negative) 
## and species on PC2 negative is mouse
## Q which are the genes that are species independent on pathologic hypertrophy

p.pcaloadings= pca_res[[1]]$rotation[,1:9] %>% 
  as_data_frame() %>% 
  mutate(gene= rownames(PCA$rotation))

saveRDS(p.pcaloadings, "data/pcaloadings.rds")


# check genes on pca loadings ---------------------------------------------


g.u<- p.pcaloadings%>% 
  arrange(desc(PC2))%>% 
  filter(PC1<0)%>%
  slice_max(order_by = PC2, n= 40)%>%
  pull(gene)

g.d<- p.pcaloadings%>% 
  arrange(desc(PC2))%>% 
  filter(PC1<0)%>%
  slice_min(order_by = PC2, n= 40)%>% pull(gene)
g.b <- 
  p.pcaloadings%>% 
  #arrange((PC1))%>% 
  filter(abs(PC2)<0.05)%>%
  slice_max(order_by = PC1, n= 40)%>% pull(gene)

cs <- c.df %>% 
  filter(!grepl("snR|prot", contrast_id) | grepl("F_snRNA_CM", contrast_id))%>% 
  pull(contrast_id) %>% 
  unique()

t(cm_res_scaled)[g.d, ]%>%
  as.data.frame()%>%
  rownames_to_column("gene")%>%
  pivot_longer(-gene)%>%
  ggplot(., 
         aes(x= gene, y= value))+
  geom_hline(yintercept = 0)+ 
  geom_boxplot(outlier.colour = NA)+
  geom_jitter(aes( #shape= animal, 
    color= name))+
  theme(axis.text.x = element_text(angle= 60 , hjust= 1))

get_box<- function(c.df, genes, c.ids){
  c.df %>% 
    filter(gene %in% c( genes),
           contrast_id %in% c.ids)%>% 
    mutate(animal= !grepl("hs", contrast_id))%>%
    mutate(gene= factor(gene, levels= c(genes)))  %>%
    group_by(contrast_id)%>%
    mutate(scale_lfc= scale(logFC))%>%
    ggplot(., 
           aes(x= gene, y= scale_lfc))+
    geom_hline(yintercept = 0)+ 
    geom_boxplot(outlier.colour = NA)+
    geom_jitter(aes( shape= animal, 
                     color= contrast_id))+
    theme(axis.text.x = element_text(angle= 60 , hjust= 1))
  
}

get_box(c.df , g.b, contrast_cm)
get_hmap2(c.df, g.b, contrast_cm)

get_box(c.df , g.d, contrast_cm)
get_hmap2(c.df, g.d, contrast_cm)


c.df %>% 
  mutate(animal= !grepl("hs", contrast_id))%>%
  filter(gene %in% c( g.u),
         contrast_id %in% cs)%>% 
  mutate(gene= factor(gene, levels= c(g.u)))%>%
  ggplot(., 
         aes(x= gene, y= logFC))+
  geom_boxplot(outlier.colour = NA)+
  geom_hline(yintercept = 0)+ 
  geom_jitter(aes( shape= animal, 
                   color= contrast_id))+
  theme(axis.text.x = element_text(angle= 60 , hjust= 1))

##add hmap

get_hmap2<- function(c.df, genes, c.ids){
  
  mat= c.df %>% 
    dplyr::select(contrast_id, gene,  logFC)%>%
    filter( contrast_id %in% c.ids)%>%
    filter(gene %in% genes)%>%
    #contrast_id %in% query_contrasts)%>%
    pivot_wider(names_from = contrast_id, values_from = logFC, values_fn = mean)%>%
    as.data.frame()%>%
    filter(!is.na(gene))%>%
    column_to_rownames("gene")
  col_names_plot= c(genes)
  
  na_sums_per_row <- apply(mat, 1, function(row) sum(is.na(row)))
  #mat_subset= mat[na_sums_per_row < 0.5 *  ncol(mat),]
  mat_subset<- mat
  x= rownames(mat_subset)
  x[!x %in% col_names_plot] <- ""
  dim(mat_subset)
  hmap_top <- ComplexHeatmap::Heatmap(t(mat_subset),
                                      #rect_gp = gpar(fill = "grey", lwd = 1),
                                      name = "logFC", 
                                      na_col = "black",
                                      border_gp = gpar(col = "black", lty = 1),
                                      cluster_columns = T,
                                      cluster_rows= T, 
                                      row_names_centered = TRUE,
                                      show_row_dend = T,
                                      show_column_dend = F,
                                      column_labels = x,
                                      column_names_side = "top",
                                      column_names_gp = gpar(fontsize = 8),
                                      row_names_side = "left",
                                      row_dend_side = "left"
                                      
  )
  hmap_top
}

get_hmap2(c.df, g.u, cs)
get_hmap2(c.df, g.d, cs)
#geom_hline(yintercept = 0)
p.pcaloadings%>%
  ggplot(., aes(x= PC3, y= PC4))+
  geom_point()

p.pcaloadings <-
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

# generate signatures -----------------------------------------------------

source("sub/global.R")
source("sub/helper.R")
library(ComplexUpset)
library(survcomp)
library(UpSetR)

pull_fisher_p<- function(query_contrasts, c.df, missing_prop= 90){
  
  # update the pval
  c.df2<- c.df%>% 
    group_by(contrast_id)%>%
    mutate(pval_min = min(pval[pval > 0]),  ## modify FDR values that are INF to smallest non-zero FDR within that contrast
           pval_mod= ifelse(pval== 0, pval_min, pval))
  c.df2 
  
  min_contrasts<- round(length(query_contrasts)*(missing_prop/100))
  print(paste0("min numbers of p : ", min_contrasts))
  
  x<- c.df2 %>% 
    filter(contrast_id %in% query_contrasts)%>%
    group_by(gene)%>% 
    select(gene, pval_mod)%>%
    mutate(gene_count = as.numeric(ave(gene, gene, FUN = length)))%>%
    filter(gene_count >= min_contrasts )%>%
    mutate(fisher_p=  survcomp::combine.test(pval_mod, "fisher", na.rm = T))%>%
    distinct(gene, fisher_p)%>%
    ungroup()%>% 
    mutate(fisher_p_adj= p.adjust(fisher_p))%>%
    arrange(fisher_p_adj)%>%
    mutate(fisher.rank = rank(fisher_p_adj, ties.method = "random"))
  
}

joint_contrast_df<-c.df
query_contrasts <- contrast_patho[!grepl("fetal", contrast_patho)]
query_contrasts <- query_contrasts[!grepl("HF", query_contrasts)]

length(query_contrasts)


## select the top sig genes for downstream signature analysis
p_hyp= pull_fisher_p(query_contrasts, c.df, 80)
hist(p_hyp$fisher_p)
hist(p_hyp$fisher_p_adj)


hyp.set <- p_hyp %>% 
  filter(fisher_p_adj< 0.001 ,
         fisher.rank< 500)%>% 
  pull(gene)

## get contrast IDs
#reduce df to alpha level cut off & contrast id
contrast_df_filt1= joint_contrast_df %>% 
  filter(contrast_id %in% query_contrasts)%>%
  filter(gene %in% hyp.set)

gene_counts = contrast_df_filt1 %>% 
  distinct(gene, contrast_id)%>%
  group_by(gene)%>%
  count()

p.df <- gene_counts %>% 
  ungroup()%>%
  count(n)%>%
  mutate(selected = ifelse(n >= missing_prop, "selected", "other"), 
         n= factor(n, levels =  1:length(query_contrasts)))


##### get the intersection genes by looking for gene with times the allowed NA
intersect_genes<- gene_counts %>% filter(n>= missing_prop)%>% pull(gene)%>% unique()

# assign whether genes are consistent in direction
df.msign= contrast_df_filt1 %>%
  dplyr::select(gene, contrast_id, logFC)%>% 
  filter(gene %in% intersect_genes) %>%
  #filter(gene=="MYL6")
  group_by(gene)%>% 
  summarise(m.sign = mean(sign(logFC)))%>%
  mutate(top_ = factor(ifelse(m.sign== 1, 
                              "upregulated", 
                              ifelse(m.sign ==-1 ,
                                     "downregulated", 
                                     "inconsistent")
  ),
  levels= c("upregulated", "inconsistent", "downregulated")
  )
  )

p.bar.intersect<- df.msign %>%
  count(top_)%>%
  ggplot(., aes(x = 1, y =n, fill = top_)) +
  geom_col(color="black" )+
  geom_text(aes(label = ifelse(n >= 10, n, "")), 
            position = position_stack(vjust = 0.5), # Center the labels within the bars
            color = "white") +
  scale_fill_manual(values = rev(c("darkblue", "darkgrey", "darkred")))+
  theme_cowplot()+
  theme(axis.text.x = element_blank(), 
        axis.line = element_blank(), 
        axis.ticks.x = element_blank())+
  labs(x="", y="number of genes", fill ="")

##signature
df.full <- p_hyp%>% left_join(df.msign, by= "gene")

# extract top genes
df.full <-
  df.full%>% 
  mutate(min_FDR = min(fisher_p_adj[fisher_p_adj > 0]),  ## modify FDR values that are INF to smallest non-zero FDR within that contrast
         FDR_mod= ifelse(fisher_p_adj== 0, min_FDR, fisher_p_adj))%>%
  mutate( label = ifelse(top_ %in% c("upregulated", "downregulated") & 
                          fisher.rank< 30,
                        gene, 
                        "") 
          ) 

##1. plot full

p.pvaldist_no_label<- df.full%>% 
  ggplot(., aes(x= fisher.rank,
                y= -log10(FDR_mod), 
                color= m.sign,
                label = label))+
  geom_point()+
  geom_vline(xintercept= 500,linetype=3)+
  geom_hline(yintercept = -log10(0.0001), color="black", linetype=3)+
  scale_color_gradient2(low="darkblue", mid="white", high= "darkred")+
  theme_cowplot()+
  xlim(c(0, max(df.full$fisher.rank)))+
  scale_x_continuous(
    breaks = c(1, 500, 1000, 2500, 5000, 7500), # Specify ranks for labels
    labels = c("1", "500", "1000", "2500", "5000", "7500") # Labels for the breaks
  )+
  theme(axis.text.x=element_text(angle= 45, hjust= 1))+
  #legend.position = "none")+
  labs(x= "p value ranking", y= "-log10(p_value)",
       color= "Average sign\nof regulation"
  )
p.pvaldist_no_label
p.pvaldist<- df.full%>% 
  ggplot(., aes(x= fisher.rank,
                y= -log10(FDR_mod), 
                color= m.sign,
                label = label))+
  geom_point()+
  geom_vline(xintercept= 500,linetype=3)+
  geom_hline(yintercept = -log10(0.0001), color="black", linetype=3)+
  scale_color_gradient2(low="darkblue", mid="white", high= "darkred")+
  theme_cowplot()+
  geom_text_repel(max.overlaps = 1000, 
                  size= 2.5, 
                  #color = "black",
                  segment.alpha = 0.3,
                  min.segment.length = unit(0.01, "cm"),
                  #nudge_x = 2,
                  force= 50, 
                  force_pull = 0.01
  )+
  xlim(c(0, max(df.full$fisher.rank)))+
  ylim(c(0, 320))+
  scale_x_continuous(
    breaks = c(1, 500, 1000, 2500, 5000, 7500), # Specify ranks for labels
    labels = c("1", "500", "1000", "2500", "5000", "7500") # Labels for the breaks
  )+
  theme(axis.text.x=element_text(angle= 45, hjust= 1))+
  #legend.position = "none")+
  labs(x= "p value ranking", y= "-log10(p_value)",
       color= "Average sign\nof regulation"
  )
p.pvaldist


##2. plot heatmap of top genes

genes <- df.full%>% 
  filter(abs(m.sign)> 0.9 & fisher.rank <500)%>%
  pull(gene)

hist(df.full$m.sign, breaks = 100)

mat= c.df %>% 
  dplyr::select(contrast_id, gene,  logFC)%>%
  filter( contrast_id %in% query_contrasts)%>%
  filter(gene %in% genes)%>%
  #contrast_id %in% query_contrasts)%>%
  pivot_wider(names_from = contrast_id, values_from = logFC, values_fn = mean)%>%
  as.data.frame()%>%
  filter(!is.na(gene))%>%
  column_to_rownames("gene")
col_names_plot= c(genes)

na_sums_per_row <- apply(mat, 1, function(row) sum(is.na(row)))

x= rownames(mat)
x[!x %in% col_names_plot] <- ""

hmap_top <- ComplexHeatmap::Heatmap(t(mat),
                                    #rect_gp = gpar(fill = "grey", lwd = 1),
                                    name = "logFC", 
                                    na_col = "black",
                                    border_gp = gpar(col = "black", lty = 1),
                                    cluster_columns = T,
                                    cluster_rows= T, 
                                    row_names_centered = TRUE,
                                    show_row_dend = F,
                                    show_column_dend = F,
                                    column_labels = x,
                                    column_names_side = "top",
                                    column_names_gp = gpar(fontsize = 8),
                                    row_names_side = "left",
                                    row_dend_side = "left"
                                    
)
hmap_top


## save results

df.full%>%
  select(-label, -min_FDR, -FDR_mod)%>%
  rename(regulation = top_)%>%
  write_csv("analysis/joint_contrast_analysis/CM_HYP_SIG.csv")

pdf("figures/signature_pval_rank.pdf",
    width= 6, height= 4)
p.pvaldist
dev.off()

pdf("figures/signature_pval_rank_nolabel.pdf",
    width= 8, height= 3)
p.pvaldist_no_label
dev.off()

pdf("figures/signature_consistencytop500.pdf",
    width= 2.8, height= 3)
p.bar.intersect
dev.off()

pdf("figures/signature_hmap.pdf",
    width= 14, height= 2.3)
hmap_top
dev.off()

hmap_top
p.bar.intersect

