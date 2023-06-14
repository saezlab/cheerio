
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
           query_contrasts= c("Mm_tac_ribo_2wk", "Hs_bulk_HCMvsNF"), 
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
    
    intersect_genes= names(venns[venns == length(query_contrasts)])
    
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


#' plot contrast expression
#' 
#' @param joint_contrast_df dataframe to be plotted
#' @param x_column name of the column to be plotted on x-axis
#' @param y_column name of the column to be plotted on y-axis
#' @param fg_name name of the column that will be plotted as a facet_grid variable
#' if the name remains NULL, no facet grid plotted
#' 
#' @return A single gene plot for multiple contrasts (ggplot object)

plot_logfc_gene = function(joint_contrast_df, 
                           x_column= "contrast_id", 
                           y_column= "logFC", 
                           fg_name= NULL, 
                           gene){
  
  p1= joint_contrast_df %>%
  ggplot(aes(x = !!rlang::ensym(x_column), y = !!rlang::ensym(y_column), fill = sig)) +
    geom_hline(yintercept = 0, color= "black")+
    geom_col(width= 0.4, color ="black") +
    theme_cowplot()+
    scale_x_discrete(drop=FALSE)+
    scale_fill_manual(values = c("TRUE" = "darkgreen",
                                 "FALSE"="orange", 
                                 drop= FALSE))+
    labs(x = "experimental group", y = "log fold change", fill = "FDR<0.05") +
    theme(panel.grid.major = element_line(color = "grey",
                                          linewidth = 0.1,
                                          linetype = 1),
          panel.border = element_rect(fill= NA, linewidth=1, color= "black"), 
          panel.grid.minor = element_blank(),
          axis.text = element_text(size= 11), 
          axis.title = element_text(size= 10)) +
    coord_flip()+
    ggtitle(gene)
  
  if(!is.null(fg_name)){p1= p1+facet_grid(rows= vars(!!rlang::ensym(fg_name)), scales = "free")}
  p1
}
