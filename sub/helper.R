#' Plot a Gene Set Enrichment Analysis (GSEA) plot for a given gene signature 
#' and gene set
#' 
#' @param signature dataframe of a gene signature. The dataframe must contain at
#'  least the gene ids in a column named \code{gene} and a column (with custom 
#'  column name) containing gene level statistics (e.g. logFC or t-values). 
#'  Additional columns will be ignored.
#' @param geneset dataframe of a gene set. The dataframe must contain at
#'  least the gene ids in a column named \code{gene}. Additional columns will be
#'  ignored. Note that both dataframes \code{signature} and \code{geneset} 
#'  must contain the same gene id type (e.g. gene symbols)
#' @param gene_level_stat unquoted string indicating the name of the gene 
#' level statistic column (that is stored in the dataframe \code{signature}). 
#' 
#' @return A GSEA-plot (ggplot object)
make_gsea_plot = function(signature, geneset, gene_level_stat) {
  
  # scaling and ranking (decreasing) gene level statistic (e.g. logFC)
  stats_df = signature %>%
    rename(stat := {{gene_level_stat}}) %>%
    # rename(stat = logFC) %>% # uncommend line to test the function
    arrange(-stat) %>% 
    transmute(gene, 
              stat = stat/max(abs(stat)),
              rank = row_number()) 
  
  # convert datafrane to named list of gene ids and scaled gene level statistic
  stats_list = stats_df %>% 
    select(gene, stat) %>%
    deframe()
  
  # extract the ranks/indices of the geneset member within the ranked signature
  indices_of_geneset_member = stats_df %>%
    inner_join(geneset, by="gene") %>% 
    arrange(rank) %>%
    pull(rank)
  
  # run gsea
  gseaRes = calcGseaStat(stats_list, selectedStats = indices_of_geneset_member, 
                         returnAllExtremes = TRUE)
  
  # extrat results from gsea analysis
  bottoms = gseaRes %>% pluck("bottoms")
  tops = gseaRes %>% pluck("tops")
  
  # highest rank = number of genes in signature
  max_rank = stats_df %>% pull(rank) %>% max()
  
  # construct coordinates of enrichment plot
  xs = rep(indices_of_geneset_member, each=2)
  ys = rbind(bottoms, tops) %>% as.vector()
  
  # build plotting table
  enrichment_df = tibble(
    rank = c(0, xs, max_rank + 1), 
    running_sum = c(0, ys, 0),
    max_top = max(tops),
    min_bottom = min(bottoms)
  )
  
  # position of the geneset member in the ranked signature
  gene_positions = tibble(rank = indices_of_geneset_member)
  
  # extraction of highest absolute enrichment score and corresponding rank
  top_es = enrichment_df %>% 
    top_n(1, abs(running_sum))
  
  # plotting of enrichment curve
  p1 = ggplot(data = enrichment_df, aes(x = rank, y = running_sum)) + 
    geom_hline(aes(yintercept = max_top), 
               colour = "#A11035", linetype = "dashed") +
    geom_hline(aes(yintercept = min_bottom), 
               colour ="#A11035", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "black") +
    geom_segment(data = top_es, 
                 mapping = aes(x = rank, y=0, xend = rank, yend = running_sum), 
                 linetype = "dashed") +
    geom_path(size=1, color ="#0098A1") +
    lims(x = c(0, max_rank + 1)) +
    labs(y = "Enrichment Score") +
    background_grid(major = "y", minor = "none", size.major = 0.4) +
    theme(axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          plot.margin=unit(c(1,1,-0.25,1), "cm"),
          legend.position = "none",
          title = element_text(size=16),
          axis.text = element_text(size=12),
          axis.title = element_text(size=14)) +
    annotate("text", x=1, y=max(tops) + 0.05, 
             label = str_c("ES:", round(gseaRes$res,2), sep = " "), 
             size=4.5, hjust=0)
  
  # plotting of gene positions among the signature
  p2 = ggplot(stats_df, aes(x=rank, y=1)) +
    geom_tile(aes(color=rank)) +
    geom_segment(gene_positions, 
                 mapping = aes(x = rank, y = 1.51, xend = rank, yend = 3.51), 
                 size = 0.5, color="black", alpha=0.5) +
    scale_color_gradientn(colours = c("#A11035", "#F6A800", "#FFED00", 
                                      "#57AB27", "#00549F", "#612158"),
                          breaks = c(1, max_rank/2, max_rank), 
                          labels = c("High", "", "Low")) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(),
      plot.margin=unit(c(-0.25,1,1,1), "cm"),
      axis.line.x = element_blank(),
      title = element_text(size=14),
      axis.text = element_text(size=12),
      axis.title = element_text(size=14),
      legend.position = "bottom",
      legend.justification = "center"
    ) +
    labs(x = "Rank", y=NULL, color="Gene-level Statistic") +
    guides(color = guide_colorbar(barwidth = 15, ticks=F, title.vjust = 0.8, title.position = "top"))
  
  # combine both plots
  p = plot_grid(p1, p2, ncol=1, align="v", axis="l", rel_heights = c(2,1))
  return(p)
}





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
           query_contrasts= c("tac_ribo_2wk", "fetal_rna_fetal1", "HCMvsNF_Fibroblast"), 
           alpha= 0.05, 
           cutoff= 15){
    
    
    contrast_df_filt= joint_contrast_df %>% 
      filter(contrast_id %in% query_contrasts)%>%
      filter(FDR< alpha)
    
    # get overview of intersect
    venns= table(contrast_df_filt$gene)
    
    x= split(contrast_df_filt$gene, contrast_df_filt$contrast_id)
    x= lapply(x, unique)
    
    p.venn= plot(euler(x, shape = "ellipse"), quantities = TRUE)
    
    p.hist= enframe(venns) %>% ggplot(., aes(x= factor(value)))+
      geom_histogram(stat="count")+
      labs(x= "number of contrasts reporting gene")
    
    
    intersect_genes= names(venns[venns == length(query_contrasts)])
    
    df.msign= contrast_df_filt %>%
      select(gene, contrast_id, logFC)%>% 
      filter(gene %in% intersect_genes) %>%
      group_by(gene)%>% 
      summarise(m.sign = mean(sign(logFC)))%>%
      mutate(top_ = ifelse(m.sign== 1, "upregulated", 
                           ifelse(m.sign ==-1 , "downregulated", "inconsistent")))
      
    top_up = df.msign %>%
      filter(top_== "upregulated")%>% pull(gene)
    
    top_dn = df.msign %>%
      filter(top_== "downregulated")%>% pull(gene)
    
    x= split(df.msign$gene, df.msign$top_)
    x= lapply(x, unique)  
    p.venn2= plot(euler(x, shape = "ellipse"), quantities = TRUE)
    
    p.int=  cowplot::plot_grid( p.venn,p.venn2,
                               ncol = 2, labels = c("Contrast intersection(s)", 
                                                   "Intersection breakdown"))
    
    # calculate the median normalized rank across contrasts
      df.median= joint_contrast_df %>%
        filter(contrast_id %in% query_contrasts,
               gene %in% intersect_genes)%>% 
        group_by(gene)%>%
        summarise(m.r= mean(ranks3))
      
      df.full= df.median %>%
        arrange(desc(m.r)) %>%
        left_join(contrast_df_filt%>%
                    select(gene, logFC, FDR, contrast_id), by= "gene")
      
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
        
    }else{ 
      p.top_genes= NULL
      
    }
    
      
    return(list(p.hist=p.int, 
                genes= list("i"= intersect_genes, 
                            "u"= top_up, 
                            "d"= top_dn),
                df= df.full, 
                p.top_genes= p.top_genes
    ))
  }
