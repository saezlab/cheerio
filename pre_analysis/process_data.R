
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
  mutate(tp= "fetal2", 
         modal= "rna",
         model ="fetal") 
fetal2= readRDS(file = "data/fetalDEgenes_PRJNA522417.rds")%>% 
  mutate(tp= "fetal1", 
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

saveRDS(obj.list, "data/GEX.list.hypertrophy.rds")
saveRDS(contrasts, "data/contrasts.hypertrophy.rds")




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

p1= contrasts %>%
  filter(model!= "fetal")%>%
  select(-PValue, -FDR)%>%
  pivot_wider(names_from= modal, values_from = logFC, values_fn= mean)%>%
  mutate(labels= ifelse(MgiSymbol %in% genes, MgiSymbol, ""))%>%
  ggplot(aes(x= rna, y= ribo, color = labels))+
  facet_grid(rows= vars(model), 
             cols= vars(tp))+
  geom_point(show.legend = T)+
  #ggrepel::geom_label_repel(mapping= aes(label =labels ), max.overlaps = 1000, show.legend = F)+
  theme(panel.grid.major = element_line(color = "grey",
                                          size = 0.1,
                                          linetype = 1))

p1
p2= ggplotly(p1, tooltip = c("MgiSymbol"))
toWebGL(p2)
partial_bundle(p2)




# Chaffin_ HCM DCM atlas ----------------------------------------------------------------------

sc.gex= read.csv("data/single_geneset.csv")%>% as_tibble()

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
