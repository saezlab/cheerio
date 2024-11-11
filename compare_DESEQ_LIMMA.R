
# DESEQ -------------------------------------------------------------------

library(tidyverse)

#setwd(dir = "~/R-projects/Collaborations/shiny_hypertophy/")

# mouse_hypertrophy data ----------------------------------------------------------------------
obj.list =
  list("TAC"= 
         list("RNA"= 
                list("2d"= readxl::read_xlsx("raw_data/2d Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_rna_only_tac2d_list.xlsx")%>% 
                       mutate(tp= "2d", 
                              modal= "rna",
                              model ="tac") ,
                     "2wk"= readxl::read_xlsx("raw_data/2wk Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_rna_tac2wk_list.xlsx")%>% 
                       mutate(tp= "2wk", 
                              modal= "rna",
                              model ="tac") 
                ),
              "RIBO"=
                list("2d"= readxl::read_xlsx("raw_data/2d Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_ribo_only_tac2d_list.xlsx")%>% 
                       mutate(tp= "2d", 
                              modal= "ribo",
                              model ="tac") ,
                     "2wk"= readxl::read_xlsx("raw_data/2wk Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_ribo_tac2wk_list.xlsx")%>% 
                       mutate(tp= "2wk", 
                              modal= "ribo",
                              model ="tac") 
                )
         ),
       "Swim"= 
         list("RNA"= 
                list("2d"= readxl::read_xlsx("raw_data/2d Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_rna_swim2d_sedentary_list.xlsx")%>% 
                       mutate(tp= "2d", 
                              modal= "rna",
                              model ="swim") ,
                     "2wk"= readxl::read_xlsx("raw_data/2wk Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_rna_swim2wk_sedentary_list.xlsx")%>% 
                       mutate(tp= "2wk", 
                              modal= "rna",
                              model ="swim") 
                ),
              "RIBO"=
                list("2d"= readxl::read_xlsx("raw_data/2d Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_ribo_swim2d_sedentary_list.xlsx")%>% 
                       mutate(tp= "2d", 
                              modal= "ribo",
                              model ="swim"),
                     "2wk"= readxl::read_xlsx("raw_data/2wk Listen/DEG_CPM1_non_protein_coding_transcripts_still_included_normalized_ribo_swim2wk_sedentary_list.xlsx")%>% 
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
        dplyr::select(logFC, PValue, FDR,MgiSymbol, GeneDescription, tp, modal, model)
    })%>%
      do.call(rbind,.)
  })%>%
    do.call(rbind,.)
})%>%
  do.call(rbind,.)

## LIMMA
contrasts_limmma <- readRDS("data/contrasts_mm_limma.rds")

contrasts
contrast_d = contrasts%>%
  mutate(contrast_id = paste(model, modal, tp, sep= "_"),
         gene = MgiSymbol)%>%
  dplyr::select(contrast_id, gene, logFC, FDR)%>%
  mutate(contrast_id= str_replace_all(contrast_id, "_fetal_rna", ""), 
         contrast_id = ifelse(!grepl("fetal", contrast_id), paste0("Mm_", contrast_id), contrast_id),
         contrast_id= str_replace_all(contrast_id, "fetal_rna_", ""),
         meth = "d")

unique(contrast_d$contrast_id)

contrast_l = contrasts_limmma%>%
  mutate(contrast_id = paste(model, modal, tp, sep= "_"),
         FDR = adj.P.Val)%>%
  dplyr::select(contrast_id, gene, logFC, FDR)%>%
  mutate(contrast_id= str_replace_all(contrast_id, "_fetal_rna", ""), 
         contrast_id = ifelse(!grepl("fetal", contrast_id), paste0("Mm_", contrast_id), contrast_id),
         contrast_id= str_replace_all(contrast_id, "fetal_rna_", ""),
         meth= "l")
intersect(unique(contrast_l$contrast_id) ,
          unique(contrast_d$contrast_id)
          )
x = rbind(contrast_l, contrast_d)%>%
  mutate(gene = str_to_title(gene))%>%
  group_by(gene, contrast_id)%>%
  pivot_wider(names_from = meth, values_from = logFC,values_fn = mean, -FDR)
unique(x$contrast_id)

x%>%
  ggplot(.,aes(x= l, d))+
  geom_point()+
  facet_wrap(~contrast_id)

rbind(contrast_l, contrast_d)%>%
  dplyr::group_by(contrast_id, gene, FDR, meth) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 

rbind(contrast_l, contrast_d)%>%
  mutate(gene = str_to_title(gene))%>%
  group_by(meth, contrast_id)%>%
  summarise(n= n())%>%
  ggplot(., aes(x= contrast_id, y= n, fill = meth ))+
  geom_col(position = "dodge")
