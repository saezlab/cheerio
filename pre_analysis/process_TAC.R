## ---------------------------
## Author: Jan D. Lanzer
##
## Date Created: 2023-05-26
##
## Copyright (c) Jan D. Lanzer, 2023
## Email: Jan.Lanzer@bioquant.uni-heidelberg.de
##
## Purpose of script:
##
## 
# ---------------------------------------------------------------------------------------------

library(tidyverse)
library(biomaRt)
library(limma)
library(edgeR)
library(ggrepel)


files_in_subdir <- list.files("data/TAC_and_swim_Counts/", full.names = TRUE)

# Set the path to the main folder
main_folder <- "data/TAC_and_swim_Counts/"

# Get a list of all subdirectories in the main folder
subdirectories <- list.dirs(main_folder, full.names = TRUE, recursive = FALSE)

# Create an empty list to store data frames
data_frames_list <- list()

# Iterate through the subdirectories
for (subdir in subdirectories) {
  cat("Processing folder:", subdir, "\n")
  
  # List files in the subdirectory
  files_in_subdir <- list.files(subdir, full.names = TRUE)
  
  # Read each CSV file and store it in the list
  subdirectory_name <- tools::file_path_sans_ext(basename(subdir))
  data_frames_list[[subdirectory_name]] <- lapply(files_in_subdir, function(file){
    
    file_name <- tools::file_path_sans_ext(basename(file))
    shortened_file_name <- str_extract(file_name, "^[^.]+\\.[^.]+\\.[^.]+\\.[^.]+\\.[^.]+")
    x= read.csv(file, sep = '\t')
    x$file_name <- shortened_file_name
    colnames(x) <- c("ENSMUSG", "gene", "Tname", "count", "sample")
    return(as_tibble(x))
  })
}



data_frames_list$Ribo_Seq_2d_Analysis_Swim

# Create an empty list to store count matrices
count_matrices_list <- list()

# Iterate through each data frame in data_frames_list
for (subdir_name in names(data_frames_list)) {
  cat("Processing folder:", subdir_name, "\n")
  
  # Extract the data frame from data_frames_list
  current_df <- data_frames_list[[subdir_name]]
  current_df <-current_df %>% do.call(rbind,.)
  # Create a count matrix using tidyr and dplyr
  count_matrix <- current_df %>%
    dplyr::select(gene, count, sample) %>%
    pivot_wider(names_from = sample, values_from = count, values_fn= sum)
  
  count_matrix <- count_matrix %>% column_to_rownames("gene")
  # Store the count matrix in the list
  count_matrices_list[[subdir_name]] <- count_matrix
}

count_matrices_list$RNA_Seq_2d_Analysis_Swim

# Create an empty list to store target data frames
target_data_list <- list()

# Define possible levels for factors
possible_models <- c("tac", "swim","sedentary.swim", "sham.tac")
possible_modals <- c("rnaseq", "ribo")
possible_timepoints <- c("2w", "2d", "3h")
#sample_id_prefix <- "rep-"

# Iterate through each count matrix in count_matrices_list
for (subdir_name in names(count_matrices_list)) {
  cat("Processing folder:", subdir_name, "\n")
  
  # Extract the count matrix from count_matrices_list
  count_matrix <- count_matrices_list[[subdir_name]]
  
  # Extract information from the sample names using string matching
  sample_names <- colnames(count_matrix)
  info <- str_extract_all(sample_names, 
                          paste(c(possible_models, possible_modals, possible_timepoints, "rep-(.{3})"), 
                                collapse = "|"))  # Create a target data frame
  info
  df= as.data.frame(matrix(unlist(info),ncol = length(sample_names)))
  rownames(df)<- c("model", "tp", "modal", "sample_id")
  target= as_tibble(t(as.data.frame(df)))
  target$full_sample <- sample_names
  target$model = factor(target$model, levels= c("sedentary.swim", "sham.tac","tac", "swim"))
  # Store the target data frame in the list
  target_data_list[[subdir_name]] <- target
}

target_data_list

names(count_matrices_list) == names(target_data_list)

counts <- count_matrices_list
targets <- target_data_list$Ribo_Seq_2d_Analysis_TAC

name= "RNA_Seq_2wk_Analysis_TAC"

res= map(names(count_matrices_list), function(name){
  counts <- count_matrices_list[[name]]
  targets <- target_data_list[[name]]
  
  #reorder
  counts2<-counts[,targets$full_sample]
  
  targets<- targets%>% mutate(group= ifelse(model %in% c("sham.tac", "sedentary.swim"),
                                 "ct", 
                                 "treatment"), 
                            group=factor(group,levels=c("treatment", "ct"))
                            )
  # get grouping of the variable
  f<-targets$group
  
  dge<- DGEList(counts=counts2, group= f)#, group=group)
  keep <- filterByExpr(dge, min.count	= 10, min.total.count= 15, min.prop = 0.5)
  table(keep)
  
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  dge <- calcNormFactors(dge)
  
  v <- voom(dge, plot=TRUE)
  
  design = model.matrix(~0+f)
  #print(targets)
  #print("design: ")
  #print(design)
  fit = lmFit(v$E, design)
  
  #Define contrasts
  cont.matrix = makeContrasts(cont1 = ftreatment- fct,
                              levels=design)
  #print(cont.matrix)
  
  #Empirical Bayes
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)
  
  DE_results = as.data.frame(topTable(fit2,
                                      adjust.method = "BH",
                                      number = Inf)) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(desc(abs(t))) %>%
    as_tibble()%>%
    mutate(contrast= name)
  DE_results$modal <- unique(targets$modal)
  DE_results$tp  <- targets%>% filter(group== "treatment")%>% pull(tp)%>% unique()
  DE_results$model  <- targets%>% filter(group== "treatment")%>% pull(model)%>% unique()
  return(DE_results)
})


res <- res %>% do.call(rbind, .)

res2<- res %>%
  mutate(modal = str_replace_all(modal, "seq", ""),
         tp = str_replace_all(tp, "2w", "2wk"))
res2%>%
  saveRDS("data/contrasts_mm_limma.rds")
