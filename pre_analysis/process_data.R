
setwd(dir = "~/R-projects/Collaborations/shiny_hypertophy/")
library(tidyverse)

# mouse_hypertrophy data ----------------------------------------------------------------------
obj.list =
  list("TAC"= 
         list("RNA"= 
              list("2d"= readxl::read_xlsx("data/raw_data/2d Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_rna_only_tac2d_list.xlsx")%>% 
                     mutate(tp= "2d", 
                            modal= "rna",
                            model ="tac") ,
                   "2wk"= readxl::read_xlsx("data/raw_data/2wk Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_rna_tac2wk_list.xlsx")%>% 
                     mutate(tp= "2wk", 
                            modal= "rna",
                            model ="tac") 
                   ),
            "RIBO"=
              list("2d"= readxl::read_xlsx("data/raw_data/2d Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_ribo_only_tac2d_list.xlsx")%>% 
                     mutate(tp= "2d", 
                            modal= "ribo",
                            model ="tac") ,
                   "2wk"= readxl::read_xlsx("data/raw_data/2wk Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_ribo_tac2wk_list.xlsx")%>% 
                     mutate(tp= "2wk", 
                            modal= "ribo",
                            model ="tac") 
                   )
            ),
       "Swim"= 
         list("RNA"= 
                list("2d"= readxl::read_xlsx("data/raw_data/2d Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_rna_swim2d_sedentary_list.xlsx")%>% 
                       mutate(tp= "2d", 
                              modal= "rna",
                              model ="swim") ,
                     "2wk"= readxl::read_xlsx("data/raw_data/2wk Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_rna_swim2wk_sedentary_list.xlsx")%>% 
                       mutate(tp= "2wk", 
                              modal= "rna",
                              model ="swim") 
                ),
              "RIBO"=
                list("2d"= readxl::read_xlsx("data/raw_data/2d Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_ribo_swim2d_sedentary_list.xlsx")%>% 
                       mutate(tp= "2d", 
                              modal= "ribo",
                              model ="swim"),
                     "2wk"= readxl::read_xlsx("data/raw_data/2wk Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_ribo_swim2wk_sedentary_list.xlsx")%>% 
                       mutate(tp= "2wk", 
                              modal= "ribo",
                              model ="swim")
                )
         )
     )

# update colnames:

obj.list= lapply(obj.list, function(mod){
  lapply(mod, function(dat){
    #print(colnames(dat))
    lapply(dat, function(tp){
     # print(colnames(tp))
      
      colnames(tp) = str_remove_all(colnames(tp), "_RNA")
      colnames(tp) = str_remove_all(colnames(tp), "_Ribo")
      return(tp)
    })
  })
})

##extract contrast data into single tidy format
contrasts= lapply(obj.list, function(mod){
  
  lapply(mod, function(dat){
      #print(colnames(dat))
  
    lapply(dat, function(tp){
      #print(colnames(tp))
        tp %>% 
           #filter(MgiSymbol %in% gene)%>% 
           select(logFC, PValue, FDR,MgiSymbol, GeneDescription, tp, modal, model)
      })%>%
      do.call(rbind,.)
    })%>%
    do.call(rbind,.)
})%>%
  do.call(rbind,.)


## add fetal expression data:
fetal1= readRDS(file = "data/fetalDEgenes_GSE52601.rds") %>% 
  mutate(tp= "Hs_fetal_Akat14", 
         modal= "rna",
         model ="fetal") 
fetal2= readRDS(file = "data/fetalDEgenes_PRJNA522417.rds")%>% 
  mutate(tp= "Hs_fetal_Spurell22", 
         modal= "rna",
         model ="fetal") 

fetal.df= 
rbind(fetal1, fetal2)%>% 
  mutate(gene = str_to_title(gene), ## change here for adequate homologue mapping
         GeneDescription= "")%>%
  rename(MgiSymbol= gene, 
         PValue= P.Value,
         FDR= adj.P.Val)%>% 
  select(colnames(contrasts))

contrasts= rbind(fetal.df, contrasts)


## add Pe



## save
saveRDS(obj.list, "data/GEX.list.hypertrophy.rds")
saveRDS(contrasts, "data/contrasts.hypertrophy.rds")


contrasts




# plot testing --------------------------------------------------------------------------------


## example  contrast plots 

p.volcs= map(unique(contrasts$modal), function(x){
  contrasts %>%
    filter(modal==x)%>% 
    ggplot(aes(y = -log10(PValue), x = logFC, color = MgiSymbol, shape= model)) +
    #geom_boxplot()+
    geom_point( size = 4, alpha = 0.6) +
    theme_classic() +
    labs(x = "logFC", y = "-log10(p-value)", color = "gene") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size= 11)) +
    #geom_hline(yintercept = 0, color = "grey", linetype = 2) +
    # geom_vline(xintercept = length(HFgenes_up)+0.5, color = "black", linetype =1) +
    ggtitle(x)
  
})

p.volcs
genes=c("Nppa", "Nppb")
pls= map(genes, function(x){
  contrasts %>% 
    filter(MgiSymbol==x)%>%
    mutate(exp.group= paste(modal, tp, sep= "_"),
           sig= FDR<0.05)%>%
    ggplot(aes(x = exp.group, y = logFC, fill = sig)) +
    facet_grid(rows= vars(model))+
    #geom_boxplot()+
    geom_col(width= 0.4, color ="black") +
    theme_cowplot() +
    scale_fill_manual(values = c("TRUE" = "darkgreen",
                                 "FALSE"="orange"))+
    labs(x = "experimental group", y = "log fold change", fill = "FDR<0.05") +
    theme(#panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size= 11), 
          axis.title = element_text(size= 10)) +
    coord_flip()+
    ggtitle(x)
  
})
pls  
cowplot::plot_grid(plotlist = pls)
## raw gene expression

genes= c("Nppa", "Nppb")
genes= c("NPPA", "NPPB")

obj.list$TAC$RNA$`2d`%>%
  select(-"...1"  ,
         -"logFC",
         -"logCPM",
         -"PValue",-"FDR",
         -"Ensemble_ID",
          -"GeneDescription",
         -"Biotype",
         -"tp",
         -"modal",
         -"model"  )%>% 
  pivot_longer(names_to = "sample_id", values_to= "gex", -MgiSymbol)
  


filter(MgiSymbol %in% genes)
  
colnames(obj.list$TAC$RNA$`2d`)
##

contrasts_HF

pls= map(genes, function(x){
  contrasts_HF %>% 
    filter(gene==x)%>%
    mutate( sig= adj.P.Val<0.05)%>%
    ggplot(aes(x = study , y = logFC, fill = sig)) +
    #facet_grid(rows= vars(model))+
    #geom_boxplot()+
    geom_col(width= 0.4, color ="black") +
    theme_cowplot() +
    labs(x = "experimental group", y = "log fold change", fill = "FDR<0.05") +
    theme(#panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size= 11), 
      axis.title = element_text(size= 10)) +
    coord_flip()+
    ggtitle(x)
  
})
pls  
cowplot::plot_grid(plotlist = pls)


## plot the correlation between ribo and rna
plot.df= contrasts %>%
  filter(model!= "fetal")%>%
  select(-PValue, -FDR)%>%
  pivot_wider(names_from= modal, values_from = logFC, values_fn= mean)%>%
  mutate(labels= ifelse(MgiSymbol %in% genes, MgiSymbol, "background"),
         labels= factor(labels, levels= c(genes, "background")),
         #labls= factor(labels, levels= c("", genes)),
         alphas= factor(ifelse(labels=="background", "bg","normal"))
  )%>%
  arrange(desc(labels))


  #prepare color palette:

  if(length(genes)==2){
    myColors <- c("green", "blue", "grey")
  }else if(length(genes)==1){
    myColors <- c("green", "grey")
  }else{
    myColors <- c(brewer.pal(length(genes), "Spectral"), "grey")
  }
  names(myColors) <- levels(plot.df$labels)
  
  p1= plot.df %>% 
  ggplot(aes(x= rna, y= ribo, color = labels, size= alphas, alpha= alphas))+
  facet_grid(rows= vars(model), 
             cols= vars(tp))+
  geom_hline(yintercept = 0, color= "darkgrey", size= 0.4)+
  geom_vline(xintercept = 0, color= "darkgrey", size= 0.4)+
  geom_point(show.legend = T)+
  scale_colour_manual("genes", values= myColors)+
  geom_abline(slope= 1, intercept = 0, color= "black", size= 0.4)+
  scale_alpha_manual(values=c("bg"= 0.3, "normal"= 1), guide = 'none')+
  scale_size_manual(values=c("bg"= 0.5, "normal"= 2), guide = 'none')+
  #ggrepel::geom_label_repel(mapping= aes(label =labels ), max.overlaps = 1000, show.legend = F)+
  theme(panel.grid.major = element_line(color = "grey",
                                        size = 0.1,
                                        linetype = 1))+
  labs(alpha= "")+
  xlab("logFC - transcriptome")+
  ylab("logFC - translatome")
p1




# Chaffin_ HCM DCM atlas ----------------------------------------------------------------------

sc.gex= read.csv("data/sc_gex_chaffin.csv")%>% as_tibble()

genes=c("NPPA", "NPPB")

pls= map(genes, function(x){
  sc.gex %>% 
    filter(Gene==x)%>%
    mutate(Significant= factor(ifelse(Significant==1, TRUE, FALSE)),
           Comparison = factor(Comparison, levels= c("DCMvsNF", "HCMvsNF", "DCMvsHCM")))%>%
    ggplot(aes(x = CellType, y = logFC, fill = Significant)) +
    facet_grid(rows= vars(Comparison))+
    geom_hline(yintercept = 0, color= "black")+
    geom_col(width= 0.4, color ="black") +
    theme_cowplot() +
    scale_fill_manual(values = c("TRUE" = "darkgreen",
                                  "FALSE"="orange"))+
    labs(x = "experimental group", y = "log fold change", fill = "FDR<0.05") +
    theme(#panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size= 11), 
      axis.title = element_text(size= 10)) +
    coord_flip()+
    ggtitle(x)
  
})

cowplot::plot_grid(plotlist = pls)
unique(sc.gex$Significant)


# explore PE in vitro -------------------------------------------------------------------------

df= readxl::read_excel("data/raw_data/invitro_PE/ribo_seq_nrvm_online_merged.xlsx")

df= read_csv("data/raw_data/invitro_PE/ribo_seq_nrvm_online_merged.csv")

df %>% mutate_all( as.integer)

df= df %>% as.data.frame()%>% column_to_rownames("EnsemblGenes"
)

dge<- DGEList(counts=df)#, group=group)

#detect and remove low expressed gene
keep <- filterByExpr(dge, min.count	= 5, min.total.count= 10, min.prop = 0.5)
table(keep)
?filterByExpr
dge <- dge[keep,,keep.lib.sizes=FALSE]

dge <- calcNormFactors(dge)

# use limma voom to transform count data to log2 counts per million
v <- voom(dge, plot=TRUE)
boxplot(df)

boxplot(v$E)

PCA <- prcomp(t(v$E[,]) ,center = TRUE, scale. = F)

x= cor(v$E)
pheatmap::pheatmap(x)
plot.pca = PCA$x %>%
  as.data.frame %>%
  rownames_to_column("Run") 
p.pca = ggplot(plot.pca,aes(x= PC1, y= PC2))+
  geom_point(size= 3)+
  theme_minimal()+
  labs(x= paste0("PC1 (",as.character(round(PCA$sdev[1]^2/sum(PCA$sdev^2)*100)),"%)"),
       y= paste("PC2 (",as.character(round(PCA$sdev[2]^2/sum(PCA$sdev^2)*100)),"%)"))+
  ggtitle(paste0(""))

p.pca


##

#annotate rat and mgi
library(biomaRt)


x1= read_csv("raw_data/invitro_PE/Ribo_DEG_CPM10_NRVM_online.csv")
x2= read_csv("raw_data/invitro_PE/RNA_DEG_CPM10_NRVM_online.csv")


rat_genes = unique(c(x1$Ensemble_ID, x2$Ensemble_ID_RNA))

#translate gene ID to symbol
ensembl <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl")

listFilters(ensembl)

foo <- getBM(attributes=c('ensembl_gene_id',
                          'external_gene_name'),
             #'mgi_symbol')
             filters = 'ensembl_gene_id',
             values = rat_genes,
             mart = ensembl)


x1= x1 %>%
  rename(ensembl_gene_id= Ensemble_ID)%>%
  left_join(foo)
  
x1 = x1 %>% dplyr::rename(gene= external_gene_name)%>%
  dplyr::select(gene, logFC, FDR)%>%
  mutate(#contrast_field= "Rn_invitro", 
         contrast_id = "Rn_invitro_ribo")

x2= x2 %>%
  dplyr::rename(ensembl_gene_id= Ensemble_ID_RNA)%>%
  left_join(foo)

x2 = x2 %>%
  dplyr::rename(gene= external_gene_name)%>%
  dplyr::select(gene, logFC, FDR)%>%
  mutate(#contrast_field= "Rn_invitro", 
         contrast_id = "Rn_invitro_rna")


saveRDS(object = rbind(x1, x2),file =  "raw_data/invitro_PE/contrast_df.rds")


