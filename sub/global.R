library(BiocManager)
options(repos = BiocManager::repositories())
library(rclipboard)
library(shiny)
library(shinyWidgets)
library(shinyhelper)
library(shinyjs)
library(DT)
#library(forcats)
library(scales)
library(dplyr)
library(tibble)
library(purrr)
library(ggplot2)
library(tidyr)
library(fgsea)
#library(decoupleR)
library(readr)
library(stringr)
library(shinycssloaders)
#plotting: 

library(plotly)
library(cowplot)
library(ggpmisc)
library(broom)
library(ggrepel)
library(cowplot)
library(ComplexHeatmap)
#install.packages(c("shinytest", "showimage", "webdriver"))
theme_set(theme_cowplot())


# load static data

##reheat
contrasts_HF = readRDS("app_data/study_contrasts.rds")
ranks = readRDS("app_data/study_ranks.rds")
directed_signature = readRDS("app_data/signature.rds")

#example
example_geneset = read_csv("app_data/multiple_geneset.csv", show_col_types = FALSE)

#progeny:
#prog.res= readRDS("app_data/progeny_results_all.rds")

#dorothea:
#df_tf= readRDS("app_data/dorothea_results_all.rds")

df_func= readRDS("app_data/funcomics_precalc.rds")

TFs <- unique(df_func%>% filter(database == "collectri")%>% pull(source))

#contrast query:

joint_contrast_df= readRDS( "app_data/contrasts_query_df_translated3.rds")


ipmc_data= readRDS("app_data/ipmc_data.rds")

# heart weight data:

HW_DF= readRDS("app_data/heart_weight_precalc.rds")

hcm_contrasts<- c("hs_HCMvsNF_RNA", 
                  "hs_HCMvsDCM_RNA", 
                  "hs_HCMrEFvsNF_prot",
                  "hs_HCMpEFvsNF_prot",
                  "hs_cHYPvsNF_prot"
)
contrast_ids <- joint_contrast_df$contrast_id%>% unique()


# app stuff ---------------------------------------------------------------
legend_lfc_plot= readRDS("app_data/legend_for_lfc_plot.rds")

error_text = "Queried gene(s) were not captured\nin data"

myColors1 <- c("grey", "#4E5D6C",  # Stronger blue-gray
              "#6A7F8A",  # Deep teal-gray
              "#A4827F",  # Muted terracotta
              "#3C6E71",  # Muted teal
              "#556B2F",  # Dark olive green
              "#8B687F")  # Muted mauve

myColors <- c("grey", "#4A90A4",  # Brighter teal-blue
                      "#6F9F8A",  # Light teal-green
                      "#E29F91",  # Soft coral-pink
                      "#3A6EA5",  # Rich blue
                      "#5B8C5A",  # Brighter olive green
                      "#AB82A4")  # Soft lavender
# Display the colors in R

myColors_soft <- c("#6A7D8A",  # Muted blue-gray
                   "#A0AAB2",  # Soft gray
                   "#D1C6AD",  # Soft beige
                   "#8C6E63",  # Muted brown
                   "#667C60",  # Muted olive green
                   "#A67D74")  # Muted pink-brown

myColors_full <- c(myColors_soft, myColors, myColors1)


# -------------------------------------------------------------------------
# ##testings:
# example_contrasts <- c(contrast_ids[grepl("TAC", contrast_ids)],
#                        contrast_ids[grepl("fetal", contrast_ids)],
#                        "hs_HCMvsNF_RNA", "hs_HCMvsNF_snRNA_CM")
# 
# 
# get_top_consistent_gene2(joint_contrast_df = joint_contrast_df,
#                          query_contrasts = example_contrasts, 
#                          missing_prop=8, 
#                          alpha= 0.05
#                          )
