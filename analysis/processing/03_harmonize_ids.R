# update contrast names ---------------------------------------------------

joint_df_translated<- readRDS("data/contrasts_query_df_translated2.rds")
unique(joint_df_translated$contrast_id)
#meta= read_csv("~/Downloads/contrast and data overview - Sheet1.csv")
meta <- read_csv("~/Downloads/contrast and data overview - Sheet2.csv")

## abreviate:
meta <- meta %>%
  mutate(modal = contrast_id %>% 
           str_split("_") %>%
           map_chr(~ tail(.x, 1))
  )

cell_type_dict <- list(
  "Cardiomyocyte" = "CM",
  "Fibroblast" = "FB",
  "Endothelial I" = "ECI",
  "Endothelial II" = "ECII",
  "Endothelial III" = "ECIII",
  "Pericyte" = "PC",
  "Macrophage" = "MP",
  "VSMC" = "VSMC", # Vascular Smooth Muscle Cell
  "Lymphocyte" = "Lympho",
  "Endocardial" = "Endo",
  "Adipocyte" = "Adipo",
  "Neuronal" = "Neu",
  "Lymphatic endothelial" = "LEC",
  "Mast cell" = "MC",
  "2wk" = "2w",
  "2d" = "2d",
  "HCMvsDCM" = "HCMvsDCM"
)

#  replace cell types using the dictionary
meta$modal_abbr <- 
  meta$modal %>%
  str_split("_") %>%
  map_chr(~ {
    last_part <- tail(.x, 1)
    if (last_part %in% names(cell_type_dict)) {
      cell_type_dict[[last_part]]
    } else {
      ""
    }
  })

# create new ID systematically
meta <- meta %>% mutate(contrast_id2 = paste(species, 
                                             `disease/model`,
                                             modality,
                                             #resolution,
                                             modal_abbr,
                                             sep= "_"), 
                        contrast_id2= str_remove(string = contrast_id2, 
                                                 pattern= "_$"))

meta$contrast_id

#translate the contrast query df

contrast_dic<- meta %>% distinct(contrast_id, contrast_id2)
unique(joint_df_translated$contrast_id)
unique(contrast_dic$contrast_id2)
# adding one abbreviatoiion
joint_df_translated$contrast_id<- str_replace_all(joint_df_translated$contrast_id,
                                                  pattern = "sc_", 
                                                  "singlecell_")
joint_df_translated_2 <-joint_df_translated%>% 
  left_join(contrast_dic)%>%
  ungroup()%>%
  dplyr::select(-contrast_id)%>%
  dplyr::select( contrast_id2, everything())%>%
  rename(contrast_id= contrast_id2)%>%
  filter(!is.na(contrast_id)) 

unique(joint_df_translated_2$contrast_id)
# here, we also remove ECII and ECII cluster
joint_df_translated_2<- 
  joint_df_translated_2%>% 
  filter(!grepl("ECII|ECIII", contrast_id)) # we remove these two redundatn contrasts

unique(joint_df_translated_2$contrast_id)

# now we add additional info on the mouse data
joint_df_translated_2<-
  joint_df_translated_2%>%
  mutate(
    model= factor(ifelse(grepl("TAC", contrast_id), "TAC", 
                         ifelse(grepl("swim", contrast_id ),
                                "swim",
                                ifelse(grepl("PE", contrast_id),
                                       "PE", NA))),
                  levels= c("swim", "TAC", "PE")), 
    modal = factor(ifelse(grepl("RNA", contrast_id), 
                          "transcriptome", 
                          ifelse(grepl("ribo",contrast_id), 
                                 "translatome", 
                                 "proteome") 
                          # levels=c("transcriptome",
                          #          "translateome", 
                          #          "proteome")
                          )
                   ), 
    timepoint = factor(ifelse(grepl("2d", contrast_id), "2d", 
                              ifelse(grepl("2w", contrast_id), "2w",NA)
    ), 
    levels= c("2d", "2w")
    )
  )
# 
 # joint_df_translated_2%>%
 #   distinct(modal, model, timepoint, contrast_id)%>%View()

 joint_df_translated_2 = 
   joint_df_translated_2 %>% 
   mutate(cc= ifelse(grepl("rn|mm", contrast_id), "A", "C") ,
          cc= ifelse(grepl("fetal", contrast_id), "D", cc),
          cc= ifelse(grepl("HCM", contrast_id), "B", cc),
          cc= ifelse(grepl("cHYP", contrast_id), "B", cc))%>%
   ungroup()
 
 joint_df_translated_2%>% distinct(cc, contrast_id)%>% print(n=100)
 
 joint_df_translated_2$contrast_id<-  str_replace_all(joint_df_translated_2$contrast_id, pattern = "hs_HCMvsDCM_RNA_HCMvsDCM", 
                 "hs_HCMvsDCM_RNA")
saveRDS(joint_df_translated_2, "app_data/contrasts_query_df_translated3.rds")
unique(joint_df_translated_2$contrast_id)
