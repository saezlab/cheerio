
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
get_top_consistent_gene<-
  function(joint_contrast_df,
           query_contrasts= c("Mm_TAC_RNA_2w", "hs_HCMvsNF_snRNA_CM", "hs_fetal_RNA"), 
           alpha= 0.05, 
           cutoff= 15,
           missing_prop= 0){
    
    
    contrast_df_filt= joint_contrast_df %>% 
      filter(contrast_id %in% query_contrasts)%>%
      filter(FDR< alpha)
    
    # get overview of intersect
    venns= table(contrast_df_filt$gene)
    
    x= split(contrast_df_filt$gene, contrast_df_filt$contrast_id)
    x= lapply(x, unique)
    
    p.venn= plot(euler(x, shape = "ellipse"), quantities = TRUE)
    
    intersect_genes= names(venns[venns >= (length(query_contrasts)*missing_prop/100)])
    
    df.msign= contrast_df_filt %>%
      dplyr::select(gene, contrast_id, logFC)%>% 
      filter(gene %in% intersect_genes) %>%
      group_by(gene)%>% 
      summarise(m.sign = mean(sign(logFC)))%>%
      mutate(top_ = ifelse(m.sign== 1, "upregulated", 
                           ifelse(m.sign ==-1 , "downregulated", "inconsistent")))
      
    top_up = df.msign %>%
      filter(top_== "upregulated")%>% pull(gene)
    
    top_dn = df.msign %>%
      filter(top_== "downregulated")%>% pull(gene)
    
    p.bar.intersect=
      ggplot(df.msign, aes(x= factor(top_, levels = c("upregulated", 
                                                    "downregulated", 
                                                    "inconsistent"))))+
      geom_bar()+
      labs(x= "consistency of regulation")+
      theme(axis.text.x = element_text(angle= 60, hjust= 1))
    
    p.int=  cowplot::plot_grid( p.venn,p.bar.intersect,
                               ncol = 2, labels = c("A", 
                                                   "B"), 
                               rel_widths = c(1,0.3))
    
    # calculate the median normalized rank across contrasts
      df.median= joint_contrast_df %>%
        filter(contrast_id %in% query_contrasts,
               gene %in% intersect_genes)%>% 
        group_by(gene)%>%
        summarise(m.r= mean(ranks3))
      
      df.full= df.median %>%
        arrange(desc(m.r)) %>%
        left_join(contrast_df_filt%>%
                    dplyr::select(gene, logFC, FDR, contrast_id), by= "gene")
      
    if(length(intersect_genes)!= 0){
      top_dn2= df.median%>% filter(gene %in% top_dn) %>% arrange(desc(m.r))%>% slice(1:cutoff)%>% pull(gene)
      top_up2= df.median%>% filter(gene %in% top_up) %>% arrange(m.r)%>% slice(1:cutoff)%>% pull(gene)
      
      p.top_gene_up=
        df.full %>% 
          filter(gene %in% c( top_up2))%>% 
          mutate(gene= factor(gene, levels= c(top_up2)))%>%
            ggplot(., 
                   aes(x= gene, y= logFC))+
            geom_boxplot(outlier.colour = NA)+
            geom_jitter(aes( color= contrast_id))+
            theme(axis.text.x = element_text(angle= 60 , hjust= 1))
            #geom_hline(yintercept = 0)
            
        p.top_gene_dn=
              df.full %>% 
              filter(gene %in% c(top_dn2 ))%>% 
              mutate(gene= factor(gene, levels= c(top_dn2 )))%>%
            ggplot(., 
                   aes(x= gene, y= logFC))+
              geom_boxplot(outlier.colour = NA)+
              geom_jitter(aes( color= contrast_id))+
              theme(axis.text.x = element_text(angle= 60 , hjust= 1))
              #geom_hline(yintercept = 0)      
            
      p.top_genes = cowplot::plot_grid(p.top_gene_up, p.top_gene_dn, 
                                       ncol = 1,
                                       labels= c("top upregulated", 
                                                 "top downregulated"))
      
      ##add hmap
      plot.genes= unique(c(top_dn, top_up))
     
      mat= joint_contrast_df %>% 
        dplyr::select(contrast_id, gene,  logFC)%>%
        filter(gene %in% plot.genes,
               contrast_id %in% query_contrasts)%>%
        pivot_wider(names_from = contrast_id, values_from = logFC, values_fn = mean)%>%
        as.data.frame()%>%
        filter(!is.na(gene))%>%
        column_to_rownames("gene")
      col_names_plot= c(top_dn2, top_up2)
      
      na_sums_per_row <- apply(mat, 1, function(row) sum(is.na(row)))
      mat_subset= mat[na_sums_per_row < 0.5 *  ncol(mat),]
      
      x= rownames(mat_subset)
      x[!x %in% col_names_plot] <- ""
    
      hmap_top <- ComplexHeatmap::Heatmap(t(mat_subset),
                                          #rect_gp = gpar(fill = "grey", lwd = 1),
                                          name = "logFC", 
                                          na_col = "black",
                                          border_gp = gpar(col = "black", lty = 1),
                                          cluster_columns = T,
                                          cluster_rows= T, 
                                          
                                          column_labels = x,
                                          column_names_gp = gpar(fontsize = 9),
                                          row_names_side = "left",
                                          row_dend_side = "left"
                                          
      )
        
    }else{ 
      p.top_genes= NULL
      hmap_top <- NULL
    }
    
    ## add hmap= 
    
    
    
    
    return(list(p.hist=p.int, 
                genes= list("i"= intersect_genes, 
                            "u"= top_up, 
                            "d"= top_dn),
                df= df.full, 
                p.top_genes= p.top_genes,
                hmap_top= hmap_top
    ))
  }

get_top_consistent_gene2 <-
  function(joint_contrast_df,
           query_contrasts= c("mm_TAC_RNA_2w", 
                              "hs_HCMvsNF_snRNA_CM",
                              "hs_fetal_RNA", 
                               "hs_HCMvsNF_RNA" ), 
           alpha= 0.05, 
           cutoff= 10,
           missing_prop= 1){
    
    #reduce df to alpha level cut off & contrast id
    contrast_df_filt= joint_contrast_df %>% 
      filter(contrast_id %in% query_contrasts)%>%
      filter(FDR< alpha)
    
    gene_counts = contrast_df_filt %>% 
      distinct(gene, contrast_id)%>%
      group_by(gene)%>%
        count()
    
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
        scale_fill_manual(values = c("selected" = "darkred", "other" = "darkgrey")) +
        theme_cowplot()+
        geom_text(
          aes(label = ifelse(selected=="selected", nn, ""), y = nn), 
          vjust = -0.5,  # Positioning of the label above the bar
          color = "darkred",  # Text color to match the highlighted bar
          size = 4        # Adjust text size as needed
          )+
        theme(axis.line = element_blank())+
        labs(x= "Number of contrasts", 
             y= "Number of DEGs",
             fill= "")
   p.overlaps
     ##### get the intersection genes by looking for gene with times the allowed NA
    intersect_genes<- gene_counts %>% filter(n>= missing_prop)%>% pull(gene)%>% unique()

    # assign whether genes are consistent in direction
    df.msign= contrast_df_filt %>%
      dplyr::select(gene, contrast_id, logFC)%>% 
      filter(gene %in% intersect_genes) %>%
      group_by(gene)%>% 
      summarise(m.sign = mean(sign(logFC)))%>%
      mutate(top_ = ifelse(m.sign== 1, "upregulated", 
                           ifelse(m.sign ==-1 , "downregulated", "inconsistent")))
    
     p.bar.intersect<- df.msign %>%
        count(top_)%>%
        ggplot(., aes(x = 1, y =n, fill = top_)) +
        geom_col(color="black" )+
        geom_text(aes(label = ifelse(n >= 10, n, "")), 
                  position = position_stack(vjust = 0.5), # Center the labels within the bars
                  color = "white") +
        scale_fill_manual(values = myColors)+
        theme_cowplot()+
        theme(axis.text.x = element_blank(), 
              axis.line = element_blank(), 
              axis.ticks.x = element_blank())+
        labs(x="", y="number of genes", fill ="")
    
    # combine
     p.int=  cowplot::plot_grid( p.overlaps,p.bar.intersect,
                                 ncol = 2, labels = c("A", 
                                                      "B"), 
                                rel_widths = c(2,1))
    
    # calculate the median normalized rank across contrasts,
    # this will serve as new ranking
      
    df.median= contrast_df_filt %>%
      filter(contrast_id %in% query_contrasts,
             gene %in% intersect_genes)%>% 
      group_by(gene)%>%
      summarise(m.r= mean(gene_rank_norm))
   
    #join back to the original data frame
    df.full= df.median %>%
      arrange(desc(m.r)) %>%
      left_join(contrast_df_filt%>%
                  dplyr::select(gene, logFC, FDR, contrast_id), by= "gene")%>%
      left_join(df.msign)
    
    if(length(intersect_genes)!= 0){
      
      # we prioritize those genes that are covered by all contrasts
      
      top_dn= df.full%>% 
        filter(gene %in% intersect_genes)%>%
        filter(top_== "downregulated")%>%
        arrange(desc(m.r))%>% 
        distinct(gene, m.r)%>%
        pull(gene)
      
      top_up= df.full%>% 
        filter(gene %in% intersect_genes)%>%
        filter(top_== "upregulated")%>%
        arrange((m.r))%>% 
        distinct(gene, m.r)%>%
        pull(gene)
      
      p <- lapply(list(top_up, top_dn), function(genes){
          
        df.full %>% 
          filter(gene %in% c( genes[1:cutoff]))%>% 
          mutate(gene= factor(gene, levels= c(genes)))%>%
          ggplot(., 
                 aes(x= gene, y= logFC))+
          geom_boxplot(outlier.colour = NA, width = 0.4)+
          geom_jitter(aes( color= contrast_id), width= 0.4)+
          scale_color_manual(values= myColors_soft)+
          theme(axis.text.x = element_text(angle= 60 , hjust= 1))+
          labs(x= "", color= "Contrast ID")
        
        
      })
      
      p.top_genes = cowplot::plot_grid(plotlist = p,
                                       ncol = 1,
                                       labels= c("Top upregulated", 
                                                 "Top downregulated"))
      
      ##add hmap
      plot.genes= unique(c(top_dn, top_up))
      
      # now we select all genes that are not inconsistent
      mat= df.full %>% 
        dplyr::select(contrast_id, gene,  logFC, top_)%>%
        filter(gene %in% plot.genes)%>%
        select(-top_)%>%
        pivot_wider(names_from = contrast_id, values_from = logFC, values_fn = mean)%>%
        as.data.frame()%>%
        #filter(!is.na(gene))%>%
        column_to_rownames("gene")%>% 
        as.matrix()
        
        # sums_per_row <-apply(mat, 1, function(row) sum(!is.na(row)))
        # mat= mat[sums_per_row>=2 ,]
        col_names_plot= c(top_dn[1:cutoff], top_up[1:cutoff])
        dim(mat)
      labeled_genes= rownames(mat)
      labeled_genes[!labeled_genes %in% col_names_plot] <- ""
      
      length(labeled_genes)
      
      hmap_top <- ComplexHeatmap::Heatmap(t(mat),
                                          #rect_gp = gpar(fill = "grey", lwd = 1),
                                          name = "logFC", 
                                          na_col = "black",
                                          border_gp = gpar(col = "black", lty = 1),
                                          cluster_columns = F,
                                          cluster_rows= F, 
                                          row_names_centered = TRUE,
                                          show_row_dend = F, 
                                          column_labels = labeled_genes,
                                          column_names_gp = gpar(fontsize = 7),
                                          row_names_side = "left",
                                          row_dend_side = "left"
                                          #clustering_distance_rows = "pearson",
                                          #clustering_distance_columns = "pearson"
                                          )
      hmap_top
      
    }else{ 
      error_text2 <- "There are no genes to be plotted.\nYou can select
      differnt contrasts or lower the number of required contrasts." 
      p <- ggplot()+ 
        annotate("text", x = 4, y = 25, size=6, label = error_text2)+
        theme_void()
      p.top_genes= p
      hmap_top = p
    }
    
    return(list(p.hist=p.int,
                #p.venn= p.venn,
                genes= list("i"= intersect_genes, 
                            "u"= top_up, 
                            "d"= top_dn),
                df= df.full, 
                p.top_genes= p.top_genes,
                hmap_top= hmap_top
    ))
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
                                 "FALSE"="darkgrey", 
                                 drop= FALSE))+
    labs(x = "", y = "log fold change", fill = "FDR<0.05") +
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
  
  legend <- cowplot::get_legend(p1)
  p1 <- remove_legend(p1)

  if(!is.null(fg_name)){
    
    # Set facet levels if facet_levels are provided
   
    p1= p1+facet_grid(rows= vars(!!rlang::ensym(fg_name)), scales = "free")
    if(colored_facet){
      p1 = make_colorful_facet_labels(p1, ...)
    }
    p1

  }
  
  return(list("p"= p1, "leg"= legend))
}

make_colorful_facet_labels <- function(p1,
                                        fills = c("#6457A6", "#EF767A","#FFE347","#23F0C7","yellow")
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

plot_composite_w_one_legend = function(p){
  all_p1_plots <- lapply(p, function(x) x$p)
  all_legends <- lapply(p, function(x) x$leg)
  all_legends <- all_legends[!is.na(all_legends)]
  
  # Combine plots without legends and the legend on the right
  final_plot <- cowplot::plot_grid(
    cowplot::plot_grid(plotlist = all_p1_plots), # All plots without legends
    all_legends[[1]],                                     # Single legend
    ncol = 2,                                   # Arrange side by side
    rel_widths = c(length(input$select_gene)*2, 2)         # Adjust width ratio
  )
  
  return(final_plot)
}

## this function plots a linear model with response variable 
## heart weight. Requires gex data frame
plot_hw_association_lm= function(HW_DF, 
                              genes,
                              my.formula = y~x){
  p= map(genes , function(x){
    #map(c("rna", "ribo"), function(y){
    to_plot_df<- HW_DF %>%
      filter(gene == x )
    
    if(dim(to_plot_df)[1]<1){
      error_text2 <- paste(x , "was not captured\nin data" )
      return(
        ggplot() + 
          annotate("text", x = 4, y = 25, size=8, label = error_text2) + 
          theme_void()   
      )
    }else{
      return(
      to_plot_df%>%
        mutate(exp.group =ifelse(grepl("ct", exp.group), "ct", exp.group))%>%
        ggplot(., aes(x= HW_BW, y= exp, color= model))+
        geom_point(aes(shape= exp.group), size = 3, alpha= 0.6)+
        facet_grid(rows= vars(modal), scales="free_y")+
        stat_smooth(fullrange = T, method = "lm", formula = my.formula, se = F, linewidth= 0.4) +
        # stat_poly_eq(aes(label = paste(after_stat(rr.label))), 
        #              label.x = "left", label.y = "top",
        #              formula = my.formula, parse = TRUE, size = 4)+
        stat_fit_glance(method = 'lm',
                        method.args = list(formula =my.formula),
                        geom = 'label_repel', 
                        aes(label = paste("P-val= ", 
                                          signif(..p.value.., digits = 1), sep = "")),
                        label.x = 'right',
                        label.y = "top",
                        size = 4,
                        alpha= 0.7)+
        ggtitle(x)+
        scale_color_manual(values = c("swim" = "darkblue",
                                      "tac"="darkred", 
                                      drop= FALSE))+
        labs(y= "Normalized gene expression",
             x= "Normalized heart weight", 
             shape= "Experimental\ngroup")+
        theme(panel.grid.major = element_line(color = "grey",
                                              linewidth = 0.1,
                                              linetype = 1),
              panel.border = element_rect(fill= NA, linewidth=1, color= "black"), 
              panel.grid.minor = element_blank(),
              axis.text = element_text(size= 11), 
              axis.title = element_text(size= 10)) 
        )
    }
  })
  
  p1= cowplot::plot_grid(plotlist = p)
  p1
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
          geom_point(show.legend = T,
                     size= 3)+
          ggtitle(paste0(x))+
          theme(axis.text.x= element_text(angle= 90, hjust= 1, vjust = 0.5),
                panel.grid.major = element_line(color = "gray80", size = 0.5), # Major grid lines
                panel.grid.minor = element_line(color = "gray90", size = 0.25)  # Minor grid lines
          )+
          geom_hline(yintercept = 0)+
          scale_shape_manual(values = c(15, 16))+
          #scale_size_continuous(range = c(1, 8))+
          scale_color_gradient(low= "grey", high= "red", limits = c(0, 1))+
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
