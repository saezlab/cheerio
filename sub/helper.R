
#plot full progeny matrix
plot_hmap= function(prog.matrix, 
                    title = "",
                    max.ps= NA,
                    ct_size= 11){
  plot.map= t(prog.matrix)
  
  require(circlize)
  if(is.na(max.ps)){
    col_fun = colorRamp2(c(min(plot.map), 0, max(plot.map)), c("blue", "white", "red"))
  }else{
  col_fun = colorRamp2(c(max.ps[1], 0, max.ps[2]), c("blue", "white", "red"))
  }
  
  #col_fun(seq(-3, 3))
  
  hmap= Heatmap(plot.map,cluster_rows = T,
                
                cluster_columns = F,
                name = "PROGENy \n score",
                column_names_rot = 90,
                column_title_gp = gpar(fontsize = ct_size),
                border = T,
                col = col_fun,
                  
                
                row_names_side= "left",
                #rect_gp = gpar(ol = "black", lty = 1),
                row_names_gp = gpar(fontsize = 10),
                column_names_gp = gpar(fontsize = 10),
                rect_gp = gpar(col = "darkgrey", lty = 1, size= 0.1),
                show_row_dend = FALSE, 
                column_title= title, 
                cell_fun = function(j, i, x, y, w, h, fill) {
                  if(abs(plot.map[i, j]) > 2) {
                    grid.text("*", x, y)
                  }
                })
  hmap
  
}



## contrast query function: 

get_top_consistent_gene2 <-
  function(joint_contrast_df,
           query_contrasts= c("mm_TAC_RNA_2w", 
                              "hs_HCMvsNF_snRNA_CM",
                              "hs_HCMvsNF_snRNA_Adipo",
                              "hs_fetal_RNA", 
                               "hs_HCMvsNF_RNA",
                              "hs_HCMvsNF_snRNA_Neu"), 
           alpha= 0.05, 
           cutoff= 10,
           missing_prop= 5){
    
    #reduce df to alpha level cut off & contrast id
    contrast_df_filt1= joint_contrast_df %>% 
      filter(contrast_id %in% query_contrasts)
    
    contrast_df_filt <- contrast_df_filt1 %>%
      filter(FDR< alpha)
    
    # count how many contrasts report a gene under these restrictions
    gene_counts = contrast_df_filt %>% 
      distinct(gene, contrast_id)%>%
      group_by(gene)%>%
        count()
    
    # add the cut off number (how many contrasts need to report a gene)
    p.df <- gene_counts %>% 
      ungroup()%>%
      count(n)%>%
      mutate(selected = ifelse(n >= missing_prop, "selected", "other"), 
             n= factor(n, levels =  1:length(query_contrasts)))
      
    ## PLOT #1 - barplot, showing number of contrasts
    p.overlaps = p.df %>%
        ggplot(aes(x = n, y= nn))+
        geom_col(aes(fill = selected), 
                       color = "black" # Keep the fill white to emphasize the outline
                       ) +  # Adjust the binwidth as needed
        scale_fill_manual(values = c("selected" = "#4D7298", "other" = "darkgrey")) +
        theme_cowplot()+
        geom_text(
          aes(label = ifelse(selected=="selected", nn, ""), y = nn), 
          vjust = -0.5,  # Positioning of the label above the bar
          color = "black",  # Text color to match the highlighted bar
          size = 4        # Adjust text size as needed
          )+
        theme(axis.line = element_blank())+
        labs(x= "Number of contrasts", 
             y= "Number of DEGs",
             fill= "")
   
    ## explore direction of regulation
    # get the intersection genes that are selected
    intersect_genes <- gene_counts %>% filter(n>= missing_prop)%>% pull(gene)%>% unique()

    # calculate mean sign of regulation and label clean regulated genes 
    df.msign= contrast_df_filt %>%
      dplyr::select(gene, contrast_id, logFC)%>% 
      filter(gene %in% intersect_genes) %>%
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
             )%>%
      ungroup()
    
    ## PLOT 2 - Stacked barplot to visualize the direction of regulation of the selected genes
     p.bar.intersect<- df.msign %>%
        count(top_)%>%
        ggplot(., aes(x = 1, y =n, fill = top_)) +
        geom_col(color="black" )+
        geom_text(aes(label = ifelse(n >= 10, n, "")), 
                  position = position_stack(vjust = 0.5), # Center the labels within the bars
                  color = "white") +
        scale_fill_manual(values = list("downregulated" = "darkblue", 
                                        "inconsistent"= "darkgrey", 
                                        "upregulated"="darkred"))+
        theme_cowplot()+
        theme(axis.text.x = element_blank(), 
              axis.line = element_blank(), 
              axis.ticks.x = element_blank())+
        labs(x="", y="number of genes", fill ="")
    
    # combine PLOT 1 & PLOT 2
    p.int=  cowplot::plot_grid( p.overlaps,p.bar.intersect,
                                 ncol = 2, labels = c("A", 
                                                      "B"), 
                                rel_widths = c(2,1.2))
    
    # select the consistent genes
    consistent.genes <- df.msign %>% 
      filter(top_ != "inconsistent") %>% 
      pull(gene)
    
    # procced only if there are any consistent genes
    if(length(consistent.genes)!= 0){
      
      # calculate the mean normalized rank across contrasts,
      # this will serve as new ranking
      # then join back to get lfc values per contrast (contrast_df_filt1)
      # as well as the assigned regulatory group (df.msign)
      df.full= contrast_df_filt %>%
        filter(contrast_id %in% query_contrasts,
               gene %in% consistent.genes)%>% 
        group_by(gene)%>%
        summarise(m.r= mean(gene_rank_norm))%>%
        arrange(desc(m.r)) %>%
        left_join(contrast_df_filt1%>%
                    dplyr::select(gene, logFC, FDR, contrast_id), by= "gene")%>%
        left_join(df.msign)
      
      # extract up and down regulated genes 
      top_dn= df.full%>% 
        filter(top_== "downregulated")%>%
        arrange(desc(m.r))%>% 
        distinct(gene, m.r)%>%
        pull(gene)
      
      top_up= df.full%>% 
        filter(top_== "upregulated")%>%
        arrange((m.r))%>% 
        distinct(gene, m.r)%>%
        pull(gene)
      
      p <- lapply(list(top_up, top_dn), function(genes){
        print(genes)
        if(length(genes)==0){
            # p <- ggplot()+ 
            # annotate("text", x = 4, y = 25, size=6, label = "no genes")+
            # theme_void()
            return(NULL)
        }
      
        # ensure that cutoff is not larger than number of genes  
        new_cut = min(cutoff, length(genes))
        
        #subset the data 
        # genes was orderd so the factor levels are ordered based on rank
        df.sub<- df.full %>% 
          filter(gene %in% c( genes[1:new_cut]))%>% 
          mutate(gene= factor(gene, levels= genes))
        
        # we want to ensure that the y-axis has a bit space and includes 0 
        if(mean(df.sub$logFC)>0){
          ylims= c(0,max(df.sub$logFC)+0.2)
        }else{
          ylims= c(min(df.sub$logFC)-0.2, 0)
        }
        
        ## PLOT 3 - boxplot showing logfc of top genes
        df.sub%>%
          ggplot(., 
                 aes(x= gene, y= logFC))+
          geom_boxplot(outlier.colour = NA, width = 0.4)+
          geom_jitter(aes( color= contrast_id), width= 0.2)+
          scale_color_manual(values= myColors_full)+
          theme(axis.text.x = element_text(angle= 60 , hjust= 1))+
          labs(x= "", color= "Contrast ID")+
          geom_hline(yintercept = 0, color="darkgrey")+
          ylim(ylims)
      })
     
      #extract legend to be plotted only once
      legends <- lapply(p, function(plot) {
        if (!is.null(plot)) {
          legend <- cowplot::get_legend(plot + theme(legend.position = "right"))
          return(legend)
        }
        return(NULL)
      })
      
      # Filter out NULL legends
      legends <- legends[!sapply(legends, is.null)]
      
      # join PLOT 3 together:
      p.top_genes = plot_grid(cowplot::plot_grid(remove_legend(p[[1]]), 
                                       remove_legend(p[[2]]),
                                       ncol = 1,
                                       labels= c("A", 
                                                 "B")),
                              legends[[1]], nrow= 1, rel_widths= c(1,0.2))
      
      ## Preparations for the heatmap
    
      # subset data and transform to matrix
      mat= df.full %>% 
        dplyr::select(contrast_id, gene,  logFC, top_)%>%
        filter(gene %in% consistent.genes)%>%
        select(-top_)%>%
        pivot_wider(names_from = contrast_id, values_from = logFC, values_fn = mean)%>%
        as.data.frame()%>%
        column_to_rownames("gene")%>% 
        as.matrix()%>%
        t()
      
      # this can be used to trouble shoot if NAs trouble the clustering of the 
      # heatmap:
      # sums_per_row <-apply(mat, 1, function(row) sum(!is.na(row)))
      # mat= mat[sums_per_row>=2 ,]
      
      # we create labeled_genes to only plot labels of genes
      # that fit the cutoff
      # for big heatmaps, this gets very messy so there we will plot nothing
      if(dim(mat)[2]>50){
        labeled_genes <- rep("",length(colnames(mat)))
      }else{
        labeled_genes= colnames(mat)
        #col_names_plot = c(top_dn[1:new_cut], top_up[1:new_cut])
        #labeled_genes[!labeled_genes %in% col_names_plot] <- ""
      }

    
      ## PLOT 4 - heatmap of all consistent genes with top genes labeled
      hmap_top <- ComplexHeatmap::Heatmap((mat),
                                          #rect_gp = gpar(fill = "grey", lwd = 1),
                                          name = "logFC", 
                                          na_col = "black",
                                          border_gp = gpar(col = "black", lty = 1),
                                          cluster_columns = F,
                                          cluster_rows= F, 
                                          row_names_centered = TRUE,
                                          show_row_dend = F, 
                                          column_labels = labeled_genes,
                                          column_names_gp = gpar(fontsize = 8),
                                          row_names_side = "left",
                                          row_dend_side = "left"
                                          #clustering_distance_rows = "pearson",
                                          #clustering_distance_columns = "pearson"
                                          )
      
      gene_list= list("i"= intersect_genes, 
                  "u"= top_up, 
                  "d"= top_dn)
      
      drugst_URL= generate_drugstone_url(c(top_dn[1:cutoff], top_up[1:cutoff]))
      
    }else{ 
      error_text2 <- "There are no genes to be plotted.\nYou can select
      differnt contrasts or lower the number of required contrasts." 
      p <- ggplot()+ 
        annotate("text", x = 4, y = 25, size=6, label = error_text2)+
        theme_void()
      p.top_genes= p
      hmap_top = p
      gene_list= list()
      drugst_URL=NULL
      df.full = NULL
    }
    
    return(list(p.hist=p.int,
                genes= gene_list,
                df= df.full, 
                p.top_genes= p.top_genes,
                hmap_top= hmap_top,
                drugst_URL=drugst_URL
    ))
  }


get_consistent_tfs <- function(df_func, 
                               query_contrasts= c("mm_TAC_RNA_2w", 
                                                  "hs_HCMvsNF_snRNA_CM",
                                                  "hs_fetal_RNA", 
                                                  "hs_HCMvsNF_RNA"),
                               alpha= 0.05,
                               use_FDR = F, 
                               #cutoff= 10,
                               missing_prop= 1
){
  #reduce df to alpha level cut off & contrast id
  if(use_FDR){
    column_p <- "FDR"
  }else{
    column_p <- "p_value"
  }
  
  contrast_df_filt= df_func %>% 
    filter(condition %in% query_contrasts)%>%
    filter(.data[[column_p]] < alpha, 
           database =="collectri")%>%
    group_by(source)%>%
    mutate(mean.score = median(score, na.rm  = T))
  
  gene_counts = contrast_df_filt %>% 
    distinct(source, condition)%>%
    group_by(source)%>%
    count()%>% arrange(desc(n))
  
  p.df <- gene_counts %>% 
    ungroup()%>%
    count(n)%>%
    mutate(selected = ifelse(n >= missing_prop, "selected", "other"), 
           n= factor(n, levels =  1:length(query_contrasts)))
  
  # if(max(p.df$n)< length(query_contrasts)){
  #   rbind(p.df, c(length(query_contrasts), 0))
  # }
  
  p.overlaps = p.df %>%
    ggplot(aes(x = n, y= nn))+
    geom_col(aes(fill = selected), 
             color = "black" # Keep the fill white to emphasize the outline
    ) +  # Adjust the binwidth as needed
    scale_fill_manual(values = c("selected" = "#4D7298", "other" = "darkgrey")) +
    theme_cowplot()+
    geom_text(
      aes(label = ifelse(selected=="selected", nn, ""), y = nn), 
      vjust = -0.5,  # Positioning of the label above the bar
      color = "darkred",  # Text color to match the highlighted bar
      size = 4        # Adjust text size as needed
    )+
    theme(axis.line = element_blank())+
    labs(x= "Number of contrasts", 
         y= "Number of TFs",
         fill= "")
  
  ##### get the intersection genes by looking for gene with times the allowed NA
  intersect_genes<- gene_counts %>% filter(n>= missing_prop)%>% pull(source)%>% unique()
  
  # assign whether genes are consistent in direction
  df.msign= contrast_df_filt %>%
    dplyr::select(source, condition, score)%>% 
    filter(source %in% intersect_genes) %>%
    group_by(source)%>% 
    summarise(m.sign = mean(sign(score)))%>%
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
    geom_text(aes(label = ifelse(n >= 1, n, "")), 
              position = position_stack(vjust = 0.5), # Center the labels within the bars
              color = "white") +
    scale_fill_manual(values = c("upregulated" = "darkred", 
                                 "inconsistent" = "darkgrey", 
                                 "downregulated" = "darkblue"))+
    theme_cowplot()+
    theme(axis.text.x = element_blank(), 
          axis.line = element_blank(), 
          axis.ticks.x = element_blank())+
    labs(x="", y="number of TFs", fill ="")
  
  # combine
  p.int=  cowplot::plot_grid( p.overlaps,p.bar.intersect,
                              ncol = 2, labels = c("A", 
                                                   "B"), 
                              rel_widths = c(2,1.2))
  
  tf_sets <- split(df.msign$source, df.msign$top_)
  tf_sets_cons <- c(tf_sets$upregulated, tf_sets$downregulated)
  
  if(length(tf_sets_cons)!= 0){
    
    # we prioritize those genes that are covered by all contrasts
    top_dn= contrast_df_filt%>% 
      filter(source %in% tf_sets$downregulated)%>%
      arrange((mean.score))%>%
      pull(source)%>% unique()
    
    top_up= contrast_df_filt%>% 
      filter(source %in% tf_sets$upregulated)%>%    
      arrange(desc(mean.score))%>%
      pull(source)%>% unique()
    
    p.top_genes<- contrast_df_filt%>%
      filter(source %in% tf_sets_cons)%>%
      ggplot(., 
             aes(x= reorder(source,mean.score),  y= score))+
      geom_boxplot(outlier.colour = NA, width = 0.4)+
      geom_jitter(aes( color= condition), width= 0.2, size = 4)+
      scale_color_manual(values= myColors_full)+
      theme(axis.text.x = element_text(angle= 60 , hjust= 1))+
      labs(x= "", color= "Contrast ID")+
      geom_hline(yintercept = 0, color="darkgrey")
    
    ##add hmap
    plot.genes= unique(c(top_dn, top_up))
    
    # now we select all genes that are not inconsistent
    mat= contrast_df_filt %>% 
      dplyr::select(condition, source,  score)%>%
      filter(source %in% tf_sets_cons)%>%
      pivot_wider(names_from = condition, values_from = score, values_fn = mean)%>%
      as.data.frame()%>%
      #filter(!is.na(gene))%>%
      column_to_rownames("source")%>% 
      as.matrix()
    
    mat<- mat[match( plot.genes, rownames(mat)),]

    hmap_top <- ComplexHeatmap::Heatmap(t(mat),
                                        #rect_gp = gpar(fill = "grey", lwd = 1),
                                        name = "TF activity", 
                                        na_col = "black",
                                        border_gp = gpar(col = "black", lty = 1),
                                        cluster_columns = F,
                                        cluster_rows= F, 
                                        row_names_centered = TRUE,
                                        show_row_dend = F, 
                                        #column_labels = labeled_genes,
                                        column_names_gp = gpar(fontsize = 8),
                                        row_names_side = "left",
                                        row_dend_side = "left"
                                        #clustering_distance_rows = "pearson",
                                        #clustering_distance_columns = "pearson"
    )
    
    gene_list= list("i"= intersect_genes, 
                    "u"= top_up, 
                    "d"= top_dn)
    
    
  }else{ 
    error_text2 <- "There are no genes to be plotted.\nYou can select
        differnt contrasts or lower the number of required contrasts." 
    p <- ggplot()+ 
      annotate("text", x = 4, y = 25, size=6, label = error_text2)+
      theme_void()
    p.top_genes= p
    hmap_top = p
    gene_list= list()
    drugst_URL=NULL
  }
  
  return(list(p.hist=p.int,
              tfs= tf_sets,
              df= contrast_df_filt, 
              p.top_genes= p.top_genes,
              hmap_top= hmap_top
  )
  )
}


#' plot contrast expression
#' 
#' @param joint_contrast_df dataframe to be plotted
#' @param x_column name of the column to be plotted on x-axis
#' @param y_column name of the column to be plotted on y-axis
#' @param fg_name name of the column that will be plotted as a facet_grid variable
#' if the name remains NULL, no facet grid plotted
#' 
#' @return A single gene plot for multiple contrasts (ggplot object)

plot_logfc_gene = function(red_contrast_df, 
                           x_column= "contrast_id", 
                           y_column= "logFC", 
                           y_column_label= "log fold change",
                           fg_name= NULL, 
                           max_fc, 
                           min_fc,
                           gene,
                           colored_facet =T,
                           ...){
  if(min_fc>0){min_fc= 0}
  if(max_fc<0){max_fc= 0}
  
  p1= red_contrast_df %>%
    ggplot(aes(x = !!rlang::ensym(x_column), y = !!rlang::ensym(y_column), fill = sig)) +
      geom_hline(yintercept = 0, color= "black")+
      geom_col(width= 0.4, color ="black") +
      theme_cowplot()+
      scale_x_discrete(drop=FALSE)+
      scale_fill_manual(values = c("TRUE" = "#4D7298",
                                   "FALSE"="grey", 
                                   drop= FALSE))+
      labs(x = "", y = y_column_label, fill = "FDR<0.05") +
      theme(panel.grid.major = element_line(color = "grey",
                                            linewidth = 0.1,
                                            linetype = 1),
            panel.border = element_rect(fill= NA, linewidth=1, color= "black"), 
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.text = element_text(size= 11), 
            axis.title = element_text(size= 10)) +
      coord_flip()+
      ggtitle(gene)+
      ylim(c(min_fc, max_fc))
  
  #legend <- cowplot::get_legend(p1)
  p1 <- remove_legend(p1)

  if(!is.null(fg_name)){
    
    # Set facet levels if facet_levels are provided
   
    p1= p1+facet_grid(rows= vars(!!rlang::ensym(fg_name)), scales = "free")
    if(colored_facet){
      p1 = make_colorful_facet_labels(p1, ...)
    }
    p1

  
  }
  return(p1)
}

# all_p1_plots <- lapply(p, function(x) x$p)
# all_legends <- lapply(p, function(x) x$leg)
# all_legends <- all_legends[!is.na(all_legends)]
# saveRDS(all_legends[[1]], "app_data/legend_for_lfc_plot.rds")

make_colorful_facet_labels <- function(p1,
                                        fills = c("#47E5BC", "#EF767A","#FFE347","#23F0C7","yellow")
                                        ){
  require(ggpubr)
  # this code is to make every facet label in a different color: 
  # credits to https://github.com/tidyverse/ggplot2/issues/2096
  g <- ggplot_gtable(ggplot_build(p1))
  stripr <- which(grepl('strip-r', g$layout$name))
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  #grid.draw(g)
  p1<- ggpubr::as_ggplot(g)
  return(p1)
}


## this function plots results of a precalculated linear model
## with response variable heart weight. 
plot_hw_association= function(HW_DF, 
                              genes){
  
  coef_vec<- HW_DF%>%
    filter(gene %in% genes )%>% pull(logcpm_coef)
  p= map(genes , function(x){
    #map(c("rna", "ribo"), function(y){
    
    to_plot_df<- HW_DF %>%
      filter(gene == x )
    
    if(dim(to_plot_df)[1]<1){
      error_text2 <- paste(x , "\nwas not captured\nin data" )
      return(
        ggplot() + 
          annotate("text", x = 4, y = 25, size=8, label = error_text2) + 
          theme_void()   
      )
    }else{
      return(
        g<- to_plot_df %>%
          mutate(group= paste(model, tp,  sep = "_"), 
                 significant= ifelse(FDR<0.05, "s", "ns"))%>%
          ggplot(aes(x= group, y= logcpm_coef, 
                     #size= -log10(logcpm_p_value),
                     color= R2, 
                     shape= significant))+
          geom_col(width= 0.05, color= NA, fill ="black")+
          geom_point(show.legend = T,
                     size= 3)+
          ggtitle(paste0(x))+
          theme(axis.text.x= element_text(angle= 90, hjust= 1, vjust = 0.5),
                panel.grid.major = element_line(color = "gray80", size = 0.5), # Major grid lines
                panel.grid.minor = element_line(color = "gray90", size = 0.25)  # Minor grid lines
          )+
          geom_hline(yintercept = 0)+
          scale_shape_manual(values = c(16, 17))+
          #scale_size_continuous(range = c(1, 8))+
          scale_color_gradient2(low= "darkgrey", mid= "blue", high= "red", limits = c(0, 1))+
          labs(x= "", y= "Coefficient",
               shape = "", 
               color = "R²")+
          ylim(c(min(coef_vec), max(coef_vec)))
      )
    }
  })
  
  # Extract legend from the first valid plot (that contains a legend)
  legend_plot <- p[[which(sapply(p, function(plot) inherits(plot, "ggplot")))[1]]]
  legend <- cowplot::get_legend(legend_plot)
 
  # Remove legends from individual plots
  p_no_legend <- lapply(p, remove_legend)
  
  # Combine plots without legends and the legend on the right
  final_plot <- cowplot::plot_grid(
    cowplot::plot_grid(plotlist = p_no_legend), # All plots without legends
    legend,                                     # Single legend
    ncol = 2,                                   # Arrange side by side
    rel_widths = c(length(genes)*2, 0.8)                      # Adjust width ratio
  )
  #p1= cowplot::plot_grid(plotlist = p)
  return(final_plot)
}

plot_transcipt_translat_corr= function(to_plot_df){
  to_plot_df %>% 
    #filter(model!= "PE")%>%
    arrange((labels))%>%
    drop_na()%>%
    ggplot(aes(x= transcriptome, y= translatome, color = labels, size= alphas, alpha= alphas))+
    { # Conditional facet_grid based on column presence
      if ("timepoint" %in% colnames(to_plot_df)) {
        facet_grid(rows = vars(model), cols = vars(timepoint))  # Facet by timepoint if it exists
      } else {
        facet_grid(rows = vars(model))  # Only facet by model if timepoint does not exist
      }
    } +
    geom_hline(yintercept = 0, color= "darkgrey", size= 0.4)+
    geom_vline(xintercept = 0, color= "darkgrey", size= 0.4)+
    geom_point(aes(colour=labels))+ 
    #geom_point(shape = 1, colour = "grey")+
    scale_color_manual("genes", values= myColors)+
    scale_alpha_manual(values=c("bg"= 0.3, "normal"= 1), guide = 'none')+
    geom_abline(slope= 1, intercept = 0, color= "black", size= 0.4, alpha= 0.8)+
    scale_size_manual(values=c("bg"= 0.5, "normal"= 2), guide = 'none')+
    #ggrepel::geom_label_repel(mapping= aes(label =labels ), max.overlaps = 1000, show.legend = F)+
    theme(panel.grid.major = element_line(color = "grey",
                                          linewidth = 0.1,
                                          linetype = 1),
          panel.border = element_rect(fill= NA, linewidth=1, color= "black"), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size= 11), 
          axis.title = element_text(size= 10)) +
    labs(alpha= "", color ="Genes")+
    xlab("logFC - transcriptome")+
    ylab("logFC - translatome")+
    guides(color = guide_legend(override.aes = list(size = 5))) 

}


# Function to remove legend from individual plots
remove_legend <- function(p) {
  p + theme(legend.position = "none")
}

generate_drugstone_url <- function(node_vec, identifier = "symbol", autofill_edges = FALSE, 
                                   interaction_ppi = "STRING", interaction_dpi = "DrugBank") {
  # Convert the node vector to a comma-separated string
  nodes <- paste(node_vec, collapse = ",")
  
  # Convert boolean options to lowercase strings
  autofill_edges_str <- tolower(as.character(autofill_edges))
  
  # Construct the base URL with parameters
  url <- paste0("https://drugst.one?",
                "nodes=", nodes,
                "&identifier=", identifier,
                "&autofillEdges=", autofill_edges_str,
                "&interactionProteinProtein=", interaction_ppi,
                "&interactionDrugProtein=", interaction_dpi)
  
  return(url)
}

make_nice_table <- function(df, color_column = NULL) {
  
  # Check if coloring is required
  if (!is.null(color_column)) {
    # Define breakpoints including 0
    brks <- quantile(df[[color_column]], probs = seq(0, 1, 0.05), na.rm = TRUE)
    brks <- unique(c(brks[brks < 0], 0, brks[brks > 0])) # Ensure 0 is included
    
    # Create a color palette
    colors <- colorRampPalette(c("#ADD8E6", "white", "#FFB6C1"))(length(brks) + 1)
  }
  
  # Create the datatable
  datatable <- df %>%
    DT::datatable(
      escape = FALSE,
      filter = "top",
      selection = list(target = 'row+column'),
      extensions = "Buttons",
      rownames = FALSE,
      options = list(
        scrollX = TRUE,
        autoWidth = TRUE,
        pageLength = 50, 
        dom = "Bfrtip",
        buttons = c("copy", "csv", "excel")
      )
    )
  
  # Apply coloring if required
  if (!is.null(color_column)) {
    datatable <- datatable %>%
      formatStyle(
        color_column,
        backgroundColor = styleInterval(brks, colors)
      )
  }
  
  return(datatable)
}

# Function to query IMPC API for selected genes

## from the data field 
# https://www.mousephenotype.org/help/programmatic-data-access/data-fields/

generate_solr_url <- function(base_url, 
                              genes, 
                              rows = 25,
                              fields = c("marker_symbol", "marker_accession_id",
                                         "mp_term_name", "zygosity", "sex", 
                                         "p_value")
                              ) {
  # Check for input validity
  if (length(genes) == 0) stop("Please provide at least one gene symbol.")
  
  # Construct the gene query
  gene_query <- paste(genes, collapse = "%20OR%20")
  
  # Construct the field query
  field_query <- paste(fields, collapse = ",")
  
  # Combine all into a full URL
  full_url <- paste0(
    base_url,
    "?q=marker_symbol:(", gene_query, ")",
    "&rows=", rows,
    "&wt=csv&fl=", field_query
  )
  
  return(full_url)
}

# process ipmc data
get_ipmc_data <- function(genes){
  
  # generate query URL
  full_url<- generate_solr_url("https://www.ebi.ac.uk/mi/impc/solr/genotype-phenotype/select",
                               genes, rows= 1000)
  
  # get table format via ipmc API
  ipmc<- read_csv(full_url)
  #process the links in html format 
  ipmc <- ipmc %>% 
    mutate(jax_link= paste0("https://www.informatics.jax.org/marker/",marker_accession_id ), 
           ipmc_link= paste0("https://www.mousephenotype.org/data/genes/", marker_accession_id ), 
           gene = paste0("<a href='", ipmc_link,"' target='_blank'>", marker_symbol,"</a>"),
           jax_link = paste0("<a href='", jax_link,"' target='_blank'>", "MGI_gene_detail","</a>"))
  
  # ordering and formatting
  ipmc %>% 
    mutate(Cardiac_Hypertrophy= factor(ifelse(mp_term_name %in% c("decreased heart weight", "increased heart weight"),
                                              mp_term_name, "-")), 
           Other_Phenotypes= factor(ifelse(Cardiac_Hypertrophy == "-", mp_term_name, "-")))%>%
    mutate(p_value = scientific(p_value))%>%
    select(gene,Cardiac_Hypertrophy, Other_Phenotypes, zygosity, sex, p_value, jax_link)%>%
    arrange(desc(Cardiac_Hypertrophy))
  
  
}



# wrappers  ---------------------------------------------------------------

get_tf_plot<- function(contrast_ids,
                       tfs, 
                       ...){
    # separate the modal time point and model vars: 
    to_plot_df= df_func %>% 
      filter(source %in% tfs,
             condition %in% contrast_ids)
    
    fc_vec = to_plot_df %>% 
      pull(score)
    
    p= map(tfs, function(x){
      to_plot_df= to_plot_df %>% 
        filter(source==x)
      
      #rows of plot df
      if(dim(to_plot_df)[1]<1){
        error_text2 <- paste(x , "\nwas not captured\nin data" )
        p <- ggplot()+ 
          annotate("text", x = 4, y = 25, size=6, label = error_text2)+
          theme_void()
        return(
          p
        )
      }else{
        p = plot_logfc_gene(red_contrast_df = to_plot_df, 
                            x_column = "condition", 
                            y_column = "score", 
                            y_column_label= "TF activity", 
                            #fg_name = "model",  
                            gene= x, 
                            max_fc = max(fc_vec), 
                            min_fc = min(fc_vec),
                            #colored_facet = F,
                            ...
        )
        return(p)
      }
      
      
    })
    
    
    # Combine plots without legends and the legend on the right
    final_plot <- cowplot::plot_grid(
      cowplot::plot_grid(plotlist = p), # All plots without legends
      legend_lfc_plot,                 # Single legend is loaded in global.R
      ncol = 2,                                   # Arrange side by side
      rel_widths = c(6, 1)         # Adjust width ratio
    )
    
    return(final_plot)
  
}

