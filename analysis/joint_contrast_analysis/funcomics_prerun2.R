library(tidyverse)
library(decoupleR)


collectri_hs<- decoupleR::get_collectri()
prog_hs<- decoupleR::get_progeny()

joint_contrast_df= readRDS("app_data/contrasts_query_df_translated3.rds")

# prepare contrast data in matrix format ------------------------------------------------------

min.size.param= 5

joint_contrast_df <- joint_contrast_df%>% 
  group_by(contrast_id)%>%
  #mutate(logFC_sc = logFC))%>%
  mutate(effect= logFC*-log10(FDR_mod))%>%
  mutate(ranked_effect= rank(effect))%>%
  mutate(scaled_effect= scale(effect),
         norm_ranked_effect = ranked_effect/max(ranked_effect))
  
joint_contrast_df %>% 
  filter(contrast_id== "hs_fetal_RNA")%>%
  ggplot(aes(x= norm_ranked_effect, y= logFC))+
  geom_point()

joint_contrast_df %>% 
  filter(contrast_id== "hs_fetal_RNA")%>%
  ggplot(aes(x= ranked_effect, y= scaled_effect))+
  geom_point()
joint_contrast_df %>% 
  filter(contrast_id== "hs_fetal_RNA")%>%
  ggplot(aes(x= effect, y= scaled_effect))+
  geom_point()

joint_contrast_df%>%
  group_by(contrast_id)%>%
  mutate(effect= scale(effect))%>%
  ggplot(aes(x= contrast_id, y= norm_ranked_effect))+
  geom_boxplot()

contrast_wide= joint_contrast_df %>%
  dplyr::select(norm_ranked_effect, contrast_id, gene)%>% 
  pivot_wider(names_from = contrast_id, values_from= norm_ranked_effect,values_fn= mean) 

c.oi <- unique(joint_contrast_df$contrast_id)
x=c.oi[21]

tfs <- map(c.oi, function(x){
  print(paste0("running enrichemnt on ", x))
  enrichme<- contrast_wide[,c("gene", x)]%>% drop_na()%>%
    as.data.frame()%>% 
    column_to_rownames("gene")%>%
    as.matrix()
  min.size.df<- collectri_hs %>% 
    filter(target %in% rownames(enrichme))%>%
    group_by(source)%>%
    count()
  if(max(min.size.df$n)< min.size.param){
    return(NULL)
  }else{
    return(decoupleR::run_ulm(enrichme,collectri_hs, minsize = min.size.param ))
  }
})%>% do.call(rbind,.)


tfs2 <- tfs %>% group_by(condition)%>%
  mutate(score_sc= scale(score) ,
         FDR= p.adjust(p_value))

hist(tfs$p_value)
hist(tfs2$FDR)


# 

prog <- map(c.oi, function(x){
  print(paste0("running enrichemnt on ", x))
  enrichme<- contrast_wide[,c("gene", x)]%>% drop_na()%>%
    as.data.frame()%>% 
    column_to_rownames("gene")%>%
    as.matrix()
  min.size.df<- prog_hs %>% filter(target %in% rownames(enrichme))%>%
    group_by(source)%>%
    count()
  if(max(min.size.df$n)< min.size.param){
    return(NULL)
  }else{
    return(decoupleR::run_ulm(enrichme,prog_hs, minsize = min.size.param ))
  }
})%>% do.call(rbind,.)

prog2<- prog%>%
  group_by(condition)%>%
  mutate(score_sc= scale(score) ,
         FDR= p.adjust(p_value))
prog2 %>% 
  ggplot(aes(x= condition, y= score_sc))+
  geom_boxplot()+
  coord_flip()

func.df <- rbind(prog2 %>% mutate(database= "progeny"),
                 tfs2 %>% mutate(database= "collectri"))%>% 
  mutate(sig= FDR<0.1)

func.df <- func.df %>% 
  ungroup()%>% 
  mutate(model= ifelse(grepl("TAC", condition), "TAC",
                       ifelse(grepl("swim", condition ),"swim",
                              "PE")),
         model = ifelse(grepl("snRNA", condition), 
                        sapply(str_split(condition, "_"), `[`, 2), 
                        model),
         
         model = ifelse(grepl("hs", condition) & grepl("prot|_RNA", condition), 
                        sapply(str_split(condition, "_"), `[`, 3), 
                        model)
  )%>%
  mutate(CellType= ifelse(grepl("snRNA", condition), 
                          sapply(str_split(condition, "_"), `[`, 4),
                          NA)
  )
func.df%>% distinct(condition, CellType)%>% print(n=100)

df_func = 
  func.df %>% 
  mutate(cc= ifelse(grepl("rn|mm", condition), "A", "C") ,
         cc= ifelse(grepl("fetal", condition), "D", cc),
         cc= ifelse(grepl("HCM", condition), "B", cc))%>%
  ungroup()
saveRDS(df_func, "app_data/funcomics_precalc.rds")
