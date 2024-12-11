library(BiocManager)
options(repos = BiocManager::repositories())
library(rclipboard)
library(shiny)
library(shinyWidgets)
library(shinyhelper)
library(shinyjs)
#library(bslib)
library(shinythemes)
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
##meta table
meta_table <- read_csv("app_data/meta_data.csv") %>%
  mutate(across(-`Unique samples`, factor)) %>%
  mutate(DOI = ifelse(
    DOI == "GSE141910", 
    '<a href="https://www.med.upenn.edu/magnet/" target="_blank">GSE141910</a>',  # Specific link for GSE141910
    ifelse(
      !is.na(DOI), 
      paste0('<a href="', DOI, '" target="_blank">', DOI, '</a>'), 
      DOI
    )
  ))   # Wrap DOI values in clickable HTML links

##reheat
contrasts_HF = readRDS("app_data/study_contrasts.rds")
ranks = readRDS("app_data/study_ranks.rds")
directed_signature = readRDS("app_data/signature.rds")

#example
example_geneset = read_csv("app_data/multiple_geneset.csv", show_col_types = FALSE)

# funcomics results
df_func= readRDS("app_data/funcomics_precalc.rds")

TFs <- unique(df_func%>% filter(database == "collectri")%>% pull(source))

#contrast query (main data frame with gene level statistics from each contrast)
joint_contrast_df= readRDS( "app_data/contrasts_query_df_translated3.rds")

# heart weight data (precalculated linear model results)
HW_DF= readRDS("app_data/heart_weight_precalc.rds")



# app stuff ---------------------------------------------------------------
hcm_contrasts<- c("hs_HCMvsNF_RNA", 
                  "hs_HCMvsDCM_RNA", 
                  "hs_HCMrEFvsNF_prot",
                  "hs_HCMpEFvsNF_prot",
                  "hs_cHYPvsNF_prot"
)

contrast_ids <- joint_contrast_df$contrast_id%>% unique()

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



