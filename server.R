# SERVER
server = function(input, output, session) {
  
#### Query genes ####

  ## reset the gene input button:
 observeEvent(input$reset_input, {
    updatePickerInput(session,
                      inputId = "select_gene", selected = character(0)
    )
  })

  #input = list(select_gene =c("Nppa"))
# cardiac mouse hypertrophy show correlation between transcript and translatome:  
output$cardiac_hyper_corr = renderPlot({
  if (!is.null(input$select_gene) ){
    plot.df= joint_contrast_df %>% 
      ungroup()%>%
      filter( grepl(pattern = "Mm|Rn", contrast_id))%>%
      mutate(gene= str_to_title(gene),
            model= factor(ifelse(grepl("tac", contrast_id), "tac", 
                                  ifelse(grepl("swim", contrast_id ),
                                         "swim", 
                                         "invitro")),levels= c("swim", "tac", "invitro")), 
             modality = factor(ifelse(grepl("rna", contrast_id), "transcriptome", "translatome")), 
             timepoint = factor(ifelse(grepl("2d", contrast_id), "2d", 
                                       ifelse(grepl("2wk", contrast_id), "2wk","none" ))) )%>% 
      select( model,modality, timepoint,gene, logFC)%>%
      pivot_wider(names_from= modality, values_from = logFC, values_fn= mean)%>%
      mutate(labels= ifelse(gene %in% input$select_gene, gene, "background"),
             labels= factor(labels, levels= c(input$select_gene, "background")),
             alphas= factor(ifelse(labels=="background", "bg","normal"))
      )%>%
      arrange(desc(labels))
  

      if(length(input$select_gene)==2){
        myColors <- c("green", "blue", "grey")
      }else if(length(input$select_gene)==1){
        myColors <- c("green", "grey")
      }else{
        myColors <- c(brewer.pal(length(input$select_gene), "Spectral"), "grey")
      }
    
      myColors <- c(brewer.pal(length(input$select_gene), "Spectral"), "grey")
      names(myColors) <- levels(plot.df$labels)
      
      
      p1= plot.df %>% 
        filter(model!= "invitro")%>%
        drop_na()%>%
        ggplot(aes(x= transcriptome, y= translatome, color = labels, size= alphas, alpha= alphas))+
        facet_grid(rows= vars(model), 
                   cols= vars(timepoint))+
        geom_hline(yintercept = 0, color= "darkgrey", size= 0.4)+
        geom_vline(xintercept = 0, color= "darkgrey", size= 0.4)+
        geom_point(show.legend = T)+
        scale_colour_manual("genes", values= myColors)+
        geom_abline(slope= 1, intercept = 0, color= "black", size= 0.4)+
        scale_alpha_manual(values=c("bg"= 0.3, "normal"= 1), guide = 'none')+
        scale_size_manual(values=c("bg"= 0.5, "normal"= 2), guide = 'none')+
        #ggrepel::geom_label_repel(mapping= aes(label =labels ), max.overlaps = 1000, show.legend = F)+
        theme(panel.grid.major = element_line(color = "grey",
                                              linewidth = 0.1,
                                              linetype = 1),
              panel.border = element_rect(fill= NA, linewidth=1, color= "black"), 
              panel.grid.minor = element_blank(),
              axis.text = element_text(size= 11), 
              axis.title = element_text(size= 10)) +
        labs(alpha= "")+
        xlab("logFC - transcriptome")+
        ylab("logFC - translatome")
      
      p2= plot.df %>% 
        drop_na()%>%
        filter(model== "invitro")%>%
        ggplot(aes(x= transcriptome, y= translatome, color = labels, size= alphas, alpha= alphas))+
        facet_grid(rows= vars(model))+
        geom_hline(yintercept = 0, color= "darkgrey", size= 0.4)+
        geom_vline(xintercept = 0, color= "darkgrey", size= 0.4)+
        geom_point(show.legend = T)+
        scale_colour_manual("genes", values= myColors)+
        geom_abline(slope= 1, intercept = 0, color= "black", size= 0.4)+
        scale_alpha_manual(values=c("bg"= 0.3, "normal"= 1), guide = 'none')+
        scale_size_manual(values=c("bg"= 0.5, "normal"= 2.5), guide = 'none')+
        #ggrepel::geom_label_repel(mapping= aes(label =labels ), max.overlaps = 1000, show.legend = F)+
        theme(panel.grid.major = element_line(color = "grey",
                                              linewidth = 0.1,
                                              linetype = 1),
              panel.border = element_rect(fill= NA, linewidth=1, color= "black"), 
              panel.grid.minor = element_blank(),
              axis.text = element_text(size= 11), 
              axis.title = element_text(size= 10),
              legend.position = "none") +
        labs(alpha= "")+
        xlab("logFC - transcriptome")+
        ylab("logFC - translatome")
      
      p= plot_grid(plot_grid(p2,NULL, ncol= 1),p1, ncol = 2, 
                   rel_widths = c(1,2.5),
                   labels= "AUTO")
      return(p)
  }
})

#  cardiac hypertrophy logFCs:
output$gene_expression_plots = renderPlot({
  if (!is.null(input$select_gene) )  {
    pls= map(input$select_gene, function(x){
      to_plot_df= joint_contrast_df %>% 
        filter(gene==toupper(x),
               grepl(pattern = "Mm|Rn", contrast_id))%>%
        mutate(model= factor(ifelse(grepl("tac", contrast_id), "tac", 
                              ifelse(grepl("swim", contrast_id ),
                                     "swim", 
                                     "invitro")),levels= c("swim", "tac", "invitro")))
      
      
      plot_logfc_gene(to_plot_df, fg_name = "model",  gene= x)
      
    })
    p1= cowplot::plot_grid(plotlist =  pls)
    p1
  }
})

# Mouse heart weight: 

output$heart_weight_plot= renderPlot({
  if (!is.null(input$select_gene) )  {
    plot_hw_association(HW_DF, input$select_gene)
  }
})

#IPMC_data
output$IPMC_table= DT::renderDataTable({
 ipmc_data %>%
    filter(Gene %in% (input$select_gene))%>%
    select(Gene, Allele, Cardiac_Hypertrophy, Other_Phenotypes,median_p_val, IPMC_link)%>%
    DT::datatable(
      escape=F, filter = "top",
      selection = "none",
      extensions = "Buttons",
      rownames = F,
      options = list(scrollX = T,
                     autoWidth = T,
                     dom = "Bfrtip",
                     buttons = c("copy", "csv", "excel","pdf"))
      )
})

## reheat (human HF): 
output$HFgene_regulation_boxplot = renderPlot({
  if (!is.null(input$select_gene) ) {
    pls= map(input$select_gene, function(x){
      # prep single study df
      to_plot_df= contrasts_HF %>% 
        filter(gene== toupper(x))%>%
        mutate( sig= adj.P.Val<0.05)
      # plot
      plot_logfc_gene(to_plot_df, x_column = study,  gene= x)
        })
    p1= cowplot::plot_grid(plotlist =  pls)
    p1
  }
})

# plot for rank positions
output$rank_position = renderPlotly({
  if (!is.null(input$select_gene) ) {
    sub_ranks = ranks %>%
      filter(gene %in% toupper(input$select_gene))
    
    max_rank = max(ranks$rank)
    
    rank_plot = sub_ranks %>%
      ggplot(aes(x=rank, y=1, label = gene)) +
      geom_segment(mapping = aes(y = 0.5, xend = rank, yend = 1.5), 
                   size = 0.5, color=aachen_color("red"), alpha=0.5) +
      geom_hline(yintercept = 1, size=1.5) +
      theme_classic() +
      theme(
        panel.grid = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(),
        axis.line.x = element_blank()
      ) +
      labs(x = "Rank", y=NULL) +
      xlim(1, max_rank)
    
    ggplotly(rank_plot, tooltip = c("label", "x"))
  }
})

# distribution of mean t-values 
output$mean_t_dist = renderPlotly({
  if (!is.null(input$select_gene)) {
    # density
    sub_ranks = ranks %>%
      filter(gene %in% toupper(input$select_gene))
    dens = ranks %>%
      ggplot(aes(x=mean_t, label = gene)) +
      stat_density(geom = "line") +
      geom_rug(data = sub_ranks, color=aachen_color("red"), size=1) +
      theme_classic() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "mean t-value", y="density")
    
    ggplotly(dens, tooltip = c("label"))
  }
})

## human HCM single cell
output$HF_single = renderPlot({
  if (!is.null(input$select_gene) ) {
    pls= map(input$select_gene, function(x){
      plot_df= joint_contrast_df%>%
        filter(grepl("singlecell", contrast_id), 
               !grepl("DCMvsNF", contrast_id))%>%
        mutate(CellType= factor(str_extract(contrast_id, "(?<=_)[^_]+$")), 
               Comparison = factor(str_extract(contrast_id, "(?<=_)[^_]+(?=_[^_]+$)"),
                                   levels= c( "HCMvsNF","HCMvsDCM")
                                   ))%>%
        filter(gene==toupper(x))
      plot_logfc_gene(plot_df, x= "CellType", fg_name = "Comparison",  gene= x)  
        
    })
    p1= cowplot::plot_grid(plotlist =  pls)
    p1
    
  }
})

## human HCM bulk
output$HFgene_regulation_magnet = renderPlot({
  if (!is.null(input$select_gene)) {
  pls= map(input$select_gene, function(x){
    joint_contrast_df %>% 
      filter(gene== toupper(x), 
             grepl("bulk", contrast_id ), 
             !grepl("DCMvsNF", contrast_id))%>%
      plot_logfc_gene(joint_contrast_df= .,  gene= x)
    })
  p1= cowplot::plot_grid(plotlist =  pls)
  p1
  }
})

## fetal gene expression :
output$fetal_gene_expression_plots = renderPlot({
  if (!is.null(input$select_gene) )  {
    pls= map(input$select_gene, function(x){
      joint_contrast_df %>% 
        filter(gene==toupper(x), 
               grepl("fetal", contrast_id ))%>%
        plot_logfc_gene(joint_contrast_df= .,  gene= x)
        
    })
    p1= cowplot::plot_grid(plotlist =  pls)
    p1
  }
})


### Contrast query ------------------------------------------------------------------------------
observeEvent(input$reset_input_contrasts, {
  updatePickerInput(session,
                    inputId = "select_contrast_hs", selected = character(0)
  )
  updatePickerInput(session,
                    inputId = "select_contrast_mm", selected = character(0)
  )
  updatePickerInput(session,
                    inputId = "select_contrast_hs2", selected = character(0)
  )
})

cont_res = eventReactive(input$submit_contrast, {
  res= get_top_consistent_gene(joint_contrast_df = joint_contrast_df, 
                          query_contrasts = c(input$select_contrast_mm,
                                              input$select_contrast_hs,
                                              input$select_contrast_hs2),
                          cutoff = input$cut_off_genes,
                          alpha= as.numeric(input$select_alpha)
                          
                          )
  return(res)
})

output$cq_hist= renderPlot({
  cont_res()$p.hist
})

output$cq_top=renderPlot({
  cont_res()$p.top_genes
})

output$cq_table_up= DT::renderDataTable({
 cont_res()$df %>% 
    filter(gene %in% cont_res()$genes$u)%>%
    group_by(gene)%>%
    mutate(mean_logFC= mean(logFC),
           mean_logFC= signif(mean_logFC, 3),
           m.r= signif(m.r, 3))%>%
    rename(mean_rank= m.r)%>%
    distinct(gene, mean_rank, mean_logFC)%>%
    arrange((mean_rank))%>%
    DT::datatable(
      # extensions = 'Buttons', 
      # options = list(
      #   dom = 'Bfrtip',
      #   buttons = c('copy', 'csv', 'excel', 'pdf', 'print')
      # ))
      escape=T, filter = "top",
      selection = list(target = 'column'),
      extensions = "Buttons", 
      rownames = F,
      options = list(scrollX = T,
                    autoWidth = T,
                    dom = "Bfrtip",
                    buttons = c("copy", "csv", "excel","pdf")))
})

output$cq_table_dn= DT::renderDataTable({
  cont_res()$df %>% 
    filter(gene %in% cont_res()$genes$d)%>%
    group_by(gene)%>%
    mutate(mean_logFC= mean(logFC),
           mean_logFC= signif(mean_logFC, 3),
           m.r= signif(m.r, 3))%>%
    rename(mean_rank= m.r)%>%
    distinct(gene, mean_rank, mean_logFC)%>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

output$cq_table_full= DT::renderDataTable({
  cont_res()$df %>% 
    rename(mean_rank= m.r)%>%
    mutate(logFC= signif(logFC, 3),
           mean_rank= signif(mean_rank, 3),
           FDR = scientific(FDR))%>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

output$clipup <- renderUI({
    rclipButton(
      inputId = "clipbtn_up",
      label = "Upregulated genes",
      clipText = paste0(cont_res()$genes$u, sep="\n"),
      #clipText = enframe(cont_res()$genes$u),
      icon = icon("clipboard")
    )
})

output$clipdn <- renderUI({
  rclipButton(
    inputId = "clipbtn_dn",
    label = "Downregulated genes",
    clipText = paste0(cont_res()$genes$d, sep="\n"),
    #clipText = enframe(cont_res()$genes$u),
    icon = icon("clipboard")
  )
})


#### Functional analysis ####
# progeny
observeEvent(input$reset_input_TF, {
  updatePickerInput(session,
                    inputId = "select_tf", selected = character(0)
  )
})


# ## A. nimal
output$funcA_tf= renderPlot({
  pls= map(input$select_tf, function(x){
    plot_df = df_tf$mm%>%
      filter(source ==str_to_title(x))%>%
      mutate(condition= paste(tp, modality,sep =  "_"),
             condition= ifelse(tp=="", paste0("Rn", condition),paste0("Mm_", condition) ))
    
    plot_df%>%plot_logfc_gene(., "condition", y= "score", gene = x, fg_name = "model")+
      labs(x = "experimental group", y = "TF activity score", fill = "p<0.05") 
  })
  p1= cowplot::plot_grid(plotlist =  pls)
  p1
})

output$funcA_pw= renderPlot({
  p1 =(plot_hmap(prog.res$mmRNA)+plot_hmap(prog.res$mmRibo)+plot_hmap(prog.res$rn))
  p1
})

output$funcA_tb_pw=DT::renderDataTable({
  lapply(names(prog.res[1:3]), function(x){
    #df= lapply(names(test), function(x){
    prog.res[[x]]%>%
      as.data.frame()%>%
      rownames_to_column("contrast")%>%
      pivot_longer(-contrast, names_to= "pathway", values_to= "score")%>%
      mutate(sig= ifelse(abs(score)>2, T, F),
             score = signif(score, 3))
  })%>% do.call(rbind, .)%>%
  mutate(modal = factor(ifelse(grepl("rna", contrast), "rna", "ribo")),
         model = factor(ifelse(grepl("tac", contrast), "tac",
                               ifelse(grepl("swim", contrast), "swim", "invitro"))),
         tp = factor(ifelse(grepl("2d", contrast), "2d", "2wk")),
         pathway= factor(pathway))%>%
  dplyr::select(-contrast)%>%
  DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                extensions = "Buttons", rownames = F,
                option = list(scrollX = T,
                              autoWidth = T,
                              dom = "Bfrtip",
                              buttons = c("copy", "csv", "excel")))
  })
output$funcA_tb_tf= DT::renderDataTable({
  df_tf$mm %>%
    dplyr::select(-statistic)%>%
    mutate_if(is.character, as.factor) %>%
    mutate(score = signif(score, 3),
           p_value = scientific(p_value),
    ) %>%
    #adj_pvalue = scientific(adj_pvalue),
    #source = factor(source, levels = sort(source))) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

## B
output$funcB_tf_sc=renderPlot({
  pls= map(input$select_tf, function(x){
    df_tf$hs_sc%>%
      filter(source ==toupper(x),
             condition != "DCMvs\nNF")%>%
      plot_logfc_gene(., "celltype", y_column= "score", gene = x, fg_name = "condition")+
      labs(x = "experimental group", y = "TF activity score", fill = "p<0.05") 
    })
  p1= cowplot::plot_grid(plotlist =  pls)
  p1
})

output$funcB_tf_tb_sc= DT::renderDataTable({
  df_tf$hs_sc%>%
    dplyr::select(-statistic)%>%
    mutate_if(is.character, as.factor) %>%
    mutate(score = signif(score, 3),
           p_value = scientific(p_value),
    ) %>%
    dplyr::select(source, condition, celltype, everything())%>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

output$funcB_tf_bulk=renderPlot({
  pls= map(input$select_tf, function(x){
    df_tf$hs_magnet %>%
      filter(source ==toupper(x),
             condition != "Hs_bulk_DCMvsNF")%>%
      plot_logfc_gene(., "condition", y_column= "score", gene = x)+
      labs(x = "experimental group", y = "TF activity score", fill = "p<0.05") 
    })
  p1= cowplot::plot_grid(plotlist =  pls)
  p1
})

output$funcB_tf_tb_bulk= DT::renderDataTable({
df_tf$hs_magnet%>%
    dplyr::select(-statistic)%>%
    mutate_if(is.character, as.factor) %>%
    mutate(score = signif(score, 3),
           p_value = scientific(p_value),
    ) %>%
    #adj_pvalue = scientific(adj_pvalue),
    #source = factor(source, levels = sort(source))) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

output$funcB_pw_bulk = renderPlot({
  p1 =plot_hmap(prog.res$hsmagnet)
  p1
})

output$funcB_pw_bulk_tb=DT::renderDataTable({
  prog.res$hsmagnet%>%
    as.data.frame()%>%
    rownames_to_column("contrast")%>%
    mutate_if(is.character, as.factor) %>%
    pivot_longer(-contrast, names_to= "pathway", values_to= "score")%>%
    mutate(sig= ifelse(abs(score)>2, T, F),
           score = signif(score, 3))%>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

output$funcB_pw_sc = renderPlot({
  sc.df.prg =prog.res$hsSC
  sc.hmaps= lapply(names(sc.df.prg), function(x){
    plot_hmap(sc.df.prg[[x]], x,
              max.ps = c(min(sapply(sc.df.prg, min)),
                         max(sapply(sc.df.prg, max)))
    )
  })
  names(sc.hmaps)= names(sc.df.prg)
  p1= eval(parse(text= paste(paste0("sc.hmaps$",paste0("`", names(sc.hmaps), "`")),  collapse = " + ")))
  #p1 =plot_hmap(contrast_ID()$prog$reheat)
  p1
})

output$funcB_pw_sc_tb=DT::renderDataTable({
  df= lapply(names(prog.res$hsSC), function(x){
    prog.res$hsSC[[x]] %>%
      as.data.frame()%>%
      rownames_to_column("contrast")%>%
      pivot_longer(-contrast, names_to= "pathway", values_to= "score")%>%
      mutate(celltype= factor(x),
             pathway= factor(pathway),
             contrast= factor(contrast),
             #sig= factor(sig)
      )
  })%>% do.call(rbind, .)
  df%>%
    as.data.frame()%>%
    filter(contrast!="DCMvs.NF")%>%
    mutate_if(is.character, as.factor) %>%
    mutate(sig= factor(ifelse(abs(score)>2, T, F)),
           score = signif(score, 3))%>%
    dplyr::select(contrast, celltype, everything())%>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})
##C
output$funcC_tf=renderPlot({
  pls= map(input$select_tf, function(x){
    df_tf$hs_reheat%>%
      filter(source ==toupper(x))%>%
      plot_logfc_gene(., "study", y= "score", gene = x)+
      labs(x = "experimental group", y = "TF activity score", fill = "p<0.05") 
      
  })
  p1= cowplot::plot_grid(plotlist =  pls)
  p1
})

output$funcC_pw= renderPlot({
  p1 =plot_hmap(prog.res$hsReheat)
  p1
})

output$funcC_pw_tb=DT::renderDataTable({
  prog.res$hsReheat%>%
    as.data.frame()%>%
    rownames_to_column("contrast")%>%
    mutate_if(is.character, as.factor) %>%
    pivot_longer(-contrast, names_to= "pathway", values_to= "score")%>%
    mutate(sig= ifelse(abs(score)>2, T, F),
           score = signif(score, 3))%>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

output$funcC_tf_tb= DT::renderDataTable({
  df_tf$hs_reheat
    dplyr::select(-statistic)%>%
    mutate_if(is.character, as.factor) %>%
    mutate(score = signif(score, 3),
           p_value = scientific(p_value),
    ) %>%
    #adj_pvalue = scientific(adj_pvalue),
    #source = factor(source, levels = sort(source))) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

##D
output$funcD_tf= renderPlot({
  pls= map(input$select_tf, function(x){
    df_tf$hs_fetal%>%
      filter(source ==toupper(x))%>%
      plot_logfc_gene(., "condition", y_column= "score", gene = x)+
      labs(x = "experimental group", y = "TF activity score", fill = "p<0.05") 
    })
  p1= cowplot::plot_grid(plotlist =  pls)
  p1
})

output$funcD_tf_tb= DT::renderDataTable({
  df_tf$hs_fetal%>%
    dplyr::select(-statistic)%>%
    mutate_if(is.character, as.factor) %>%
    mutate(score = signif(score, 3),
           p_value = scientific(p_value),
    ) %>%
    #adj_pvalue = scientific(adj_pvalue),
    #source = factor(source, levels = sort(source))) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

output$funcD_pw= renderPlot({
    p1 =plot_hmap(prog.res$hsfetal)
  p1
})

output$funcD_pw_tb=DT::renderDataTable({
  prog.res$hsfetal%>%
    as.data.frame()%>%
    rownames_to_column("contrast")%>%
    mutate_if(is.character, as.factor) %>%
    pivot_longer(-contrast, names_to= "pathway", values_to= "score")%>%
    mutate(sig= ifelse(abs(score)>2, T, F),
           score = signif(score, 3))%>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})



hide("loading-content", TRUE, "fade")


#### Custom Enrichment ####
# load user input (gene sets)
gs = reactive({
  if (input$take_example_data == F) {
    shinyjs::enable("user_input")
    inFile = input$user_input
    if (is.null(inFile)){
      return(NULL)
    }
    read_csv(inFile$datapath)
  } else {
    shinyjs::disable("user_input")
    example_geneset
  }
})

signature = reactive({
  switch(input$select_contrast_enrichment,
         "A. Animal models" = joint_contrast_df %>% filter(cc == "A"),
         "B. Human HCM" =joint_contrast_df %>% filter(cc == "B"),
         "C. Human HF" = joint_contrast_df%>% filter(cc == "C"),
         "D. Fetal reprogramming"= joint_contrast_df %>% filter(cc == "D"),
         )
})
#for trouble shooting:
#sig = joint_contrast_df %>% filter(cc == "A")
#gss = example_geneset

# perform GSEA with plotting
gsea_res = eventReactive(input$submit, {
  if (ncol(gs()) == 1) {
    df= map(unique(signature()$contrast_id), function(x){
      y= signature()%>% filter(contrast_id == x)
      vect= deframe(y[,c("gene", "logFC")])
      res = fgsea(stats= vect ,
                  pathways= gs()$gene, 
                  nperm = 1000)
      as.data.frame(res)%>% mutate(cc= x)
    }) %>% do.call(rbind,. )
    
    p1= df %>% 
      ggplot(., aes(x= cc, y= NES, fill = padj<0.05)) +
      geom_col()+
      scale_fill_manual(values = c("TRUE" = "darkgreen",
                                   "FALSE"="orange"))+
      coord_flip()+
      labs(x= "contrast ID", 
           y= "normalized enrichment score (NES)")
    
  } else if (ncol(gs()) == 2) {
    list_gs = gs()%>%
      split(.$geneset) %>%
      map(pull, gene)
    df= map(unique(signature()$contrast_id), function(x){
      y= signature()%>% filter(contrast_id == x)
      vect= deframe(y[,c("gene", "logFC")])
      res = fgsea(stats= vect ,
                  pathways= list_gs, 
                  nperm = 1000)
      as.data.frame(res)%>% mutate(cc= x)
    }) %>% do.call(rbind,. )
    
    p1= df %>% 
      ggplot(., aes(x= cc, y= NES, fill = padj<0.05)) +
      geom_col()+
      facet_grid(rows= vars(pathway))+
      scale_fill_manual(values = c("TRUE" = "darkgreen",
                                   "FALSE"="orange"))+
      coord_flip()+
      labs(x= "contrast ID", 
           y= "normalized enrichment score (NES)")
  }
  
df = df %>% 
  dplyr::select(cc, everything())%>%
  dplyr::rename(geneset = pathway,
                contrast_ID = cc) %>%
    as_tibble() %>%
    select(-leadingEdge) %>%
    mutate(signature = input$signature_source)
  
  list(df = df, p = p1)
})

# gsea results as table
output$gsea_res_table = DT::renderDataTable({
  gsea_res()$df%>%
  #df%>%
    mutate(NES = signif(NES, 3),
           ES = signif(ES, 3),
           pval = scientific(pval),
           padj = scientific(padj)) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "column"),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

# gsea results as plots
output$gsea_res_plots = renderPlot({
  gsea_res()$p
})



####Download center ####

# make gene contrast data accessible:
output$mouse_hypertrophyDT = DT::renderDataTable( {
  
  joint_contrast_df %>% 
    filter( grepl(pattern = "Mm|Rn", contrast_id))%>%
    mutate(model= factor(ifelse(grepl("tac", contrast_id), "tac", 
                                ifelse(grepl("swim", contrast_id ),
                                       "swim", 
                                       "invitro")),levels= c("swim", "tac", "invitro")), 
           modality = factor(ifelse(grepl("rna", contrast_id), "transcriptome", "translatome")), 
           timepoint = factor(ifelse(grepl("2d", contrast_id), "2d", 
                                     ifelse(grepl("2wk", contrast_id), "2wk","none" ))) )%>% 
    select( model,modality, timepoint,gene, logFC, FDR)%>%
    mutate(logFC = signif(logFC,3),
           FDR = scientific(FDR)
           ) %>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
  
})


output$human_HF_bulk_summDT = DT::renderDataTable({
  
  ranks %>%
    mutate(mean_lfc = signif(mean_lfc,3),
           mean_t = signif(mean_t,3),
           fisher_pvalue = scientific(fisher_pvalue)) %>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
  })

output$human_fetalDT = DT::renderDataTable({
  joint_contrast_df%>%
    select(contrast_id, gene, logFC, FDR,sig)%>%
    filter(grepl("fetal", contrast_id))%>%
    mutate(logFC = signif(logFC,3),
           FDR = scientific(FDR),
           contrast_id = as_factor(contrast_id)) %>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
  
  
})

output$human_HCMDT = DT::renderDataTable({
  
  joint_contrast_df%>% 
    filter(grepl(pattern="HCM", contrast_id) )%>%
    mutate(modality= factor(ifelse(grepl("bulk", contrast_id),"bulk", "single_cell")))%>%
    mutate(logFC = signif(logFC,3),
         FDR = scientific(FDR),
         contrast_id = as_factor(contrast_id)) %>%
    select(modality, contrast_id, gene, logFC, FDR,sig)%>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                                 extensions = "Buttons", rownames = F,
                                 option = list(scrollX = T,
                                               autoWidth = T,
                                               dom = "Bfrtip",
                                               buttons = c("copy", "csv", "excel")))
    

})


}


