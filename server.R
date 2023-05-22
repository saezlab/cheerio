# SERVER
server = function(input, output, session) {
  
#### Query genes ####

## murine cardiac Hyp:
#  cardiac hypertrophy logFCs:
output$gene_expression_plots = renderPlot({
  if (!is.null(input$select_gene) & ("Murine_Hypertrophy" %in% input$contrasts))  {
    pls= map(input$select_gene, function(x){
      contrasts %>% 
        filter(MgiSymbol==x, model != "fetal")%>%
        mutate(exp.group= paste(modal, tp, sep= "_"),
               sig= FDR<0.05)%>%
        ggplot(aes(x = exp.group, y = logFC, fill = sig)) +
        facet_grid(rows= vars(model))+
        #geom_boxplot()+
        geom_hline(yintercept = 0, color= "black")+
        geom_col(width= 0.4, color ="black") +
        theme_cowplot() +
        scale_fill_manual(values = c("TRUE" = "darkgreen",
                                     "FALSE"="orange"))+
        labs(x = "experimental group", y = "log fold change", fill = "FDR<0.05") +
        theme(panel.grid.major = element_line(color = "grey",
                                              size = 0.1,
                                              linetype = 1),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size= 11), 
          axis.title = element_text(size= 10)) +
        coord_flip()+
        ggtitle(x)
      
    })
    p1= cowplot::plot_grid(plotlist =  pls)
    p1
  }
})

# cardiac hypertrophy show correlation::  
output$cardiac_hyper_corr = renderPlot({
  if (!is.null(input$select_gene) & ("Murine_Hypertrophy" %in% input$contrasts)){
    
    p1= contrasts %>%
      filter(model!= "fetal")%>%
      select(-PValue, -FDR)%>%
      pivot_wider(names_from= modal, values_from = logFC, values_fn= mean)%>%
      mutate(labels= ifelse(MgiSymbol %in% input$select_gene, MgiSymbol, ""),
             labls= factor(labels, levels= c("", input$select_gene)),
             alphas= factor(ifelse(labels=="", "bg","normal"))
             )%>%
      ggplot(aes(x= rna, y= ribo, color = labels,size= alphas, alpha= alphas))+
      facet_grid(rows= vars(model), 
                 cols= vars(tp))+
      geom_point(show.legend = T)+
      scale_alpha_manual(values=c("bg"= 0.3, "normal"= 1))+
      scale_size_manual(values=c("bg"= 0.5, "normal"= 2))+
      #ggrepel::geom_label_repel(mapping= aes(label =labels ), max.overlaps = 1000, show.legend = F)+
      theme(panel.grid.major = element_line(color = "grey",
                                            size = 0.1,
                                            linetype = 1))+
      labs(alpha= "")
    
    p1
  }
})

## human HF: 
output$HFgene_regulation_boxplot = renderPlot({
  if (!is.null(input$select_gene) & ("Human_HF" %in% input$contrasts)) {
    pls= map(input$select_gene, function(x){
      contrasts_HF %>% 
        filter(gene== toupper(x))%>%
        mutate( sig= adj.P.Val<0.05)%>%
        ggplot(aes(x = study , y = logFC, fill = sig)) +
        #facet_grid(rows= vars(model))+
        #geom_boxplot()+
        geom_col(width= 0.4, color ="black") +
        theme_cowplot() +
        scale_fill_manual(values = c("TRUE" = "darkgreen",
                                     "FALSE"="orange"))+
        labs(x = "HF Study", y = "log fold change", fill = "FDR<0.05") +
        theme(panel.grid.major = element_line(color = "grey",
                                              size = 0.1,
                                              linetype = 1),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size= 11), 
          axis.title = element_text(size= 10)) +
        coord_flip()+
        ggtitle(x)
      })
    p1= cowplot::plot_grid(plotlist =  pls)
    p1
    
  }
})

# plot for rank positions
output$rank_position = renderPlotly({
  if (!is.null(input$select_gene) & ("Human_HF" %in% input$contrasts)) {
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
  if (!is.null(input$select_gene) & ("Human_HF" %in% input$contrasts)) {
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

##  single cell human
output$HF_single = renderPlot({
  if (!is.null(input$select_gene) & ("Human_HF" %in% input$contrasts)) {
    pls= map(input$select_gene, function(x){
      sc.gex %>% 
        filter(Gene==toupper(x))%>%
        # mutate(Significant= factor(ifelse(Significant==1, TRUE, FALSE)),
        #        Comparison = factor(Comparison, levels= c("DCMvsNF", "HCMvsNF", "DCMvsHCM")))%>%
        ggplot(aes(x = CellType, y = logFC, fill = Significant)) +
        facet_grid(rows= vars(Comparison))+
        geom_hline(yintercept = 0, color= "black")+
        geom_col(width= 0.4, color ="black") +
        theme_cowplot() +
        scale_fill_manual(values = c("TRUE" = "darkgreen",
                                     "FALSE"="orange"))+
        labs(x = "experimental group", y = "log fold change", fill = "FDR<0.05") +
        theme(#panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size= 11), 
          axis.title = element_text(size= 10)) +
        coord_flip()+
        ggtitle(x)
      
    })
    p1= cowplot::plot_grid(plotlist =  pls)
    p1
    
  }
})

## fetal gene expression :
output$fetal_gene_expression_plots = renderPlot({
  if (!is.null(input$select_gene) & ("Fetal" %in% input$contrasts))  {
    pls= map(input$select_gene, function(x){
      contrasts %>% 
        filter(MgiSymbol==x, 
               model == "fetal")%>%
        mutate(exp.group= paste(modal, tp, sep= "_"),
               sig= FDR<0.05)%>%
        ggplot(aes(x = exp.group, y = logFC, fill = sig)) +
        #geom_boxplot()+
        geom_hline(yintercept = 0, color= "black")+
        geom_col(width= 0.4, color ="black") +
        theme_cowplot() +
        scale_fill_manual(values = c("TRUE" = "darkgreen",
                                     "FALSE"="orange"))+
        labs(x = "experimental group", y = "log fold change", fill = "FDR<0.05") +
        theme(panel.grid.major = element_line(color = "grey",
                                              size = 0.1,
                                              linetype = 1),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size= 11), 
          axis.title = element_text(size= 10)) +
        coord_flip()+
        ggtitle(x)
      
    })
    p1= cowplot::plot_grid(plotlist =  pls)
    p1
  }
})

## remove input button:
observeEvent(input$reset_input, {
  updatePickerInput(session,
                    inputId = "select_gene", selected = character(0)
                    )
  #updateTextInput(session, "mytext", value = "test")
})


# # subsetted dataframe of meta analysis
# output$summary_sub = DT::renderDataTable({
#   if (!is.null(input$select_gene)) {
#     ranks %>%
#       filter(gene %in% input$select_gene) %>%
#       mutate(mean_lfc = signif(mean_lfc,3),
#              mean_t = signif(mean_t,3),
#              fisher_pvalue = scientific(fisher_pvalue)) %>%
#       DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
#                     extensions = "Buttons", rownames = F,
#                     option = list(scrollX = T, 
#                                   autoWidth = T, 
#                                   dom = "Bfrtip",
#                                   buttons = c("copy", "csv", "excel")))
#   }
# })

#### Input data ####
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

# get which signature should be used (directed vs undirected)
signature = reactive({
  
  ## add statement of 
  switch(input$contrasts_gsea,
         "directed" = directed_signature,
         "undirected" = undirected_signature)
})

# perform GSEA with plotting
gsea_res = eventReactive(input$submit, {
  if (ncol(gs()) == 1) {
    res = fgsea(list(geneset = gs()$gene), deframe(signature()), nperm = 1000)
    p = make_gsea_plot(signature(), gs(), weight)
  } else if (ncol(gs()) == 2) {
    list_gs = gs() %>%
      split(.$geneset) %>%
      map(pull, gene)
    res = fgsea(list_gs, deframe(signature()), nperm = 1000)
    
    plot_df = gs() %>%
      nest(set = gene) %>%
      mutate(plot = pmap(., .f = function(geneset, set, ...) {
        make_gsea_plot(signature(), set, weight)
      }))
    
    p = plot_grid(plotlist = plot_df$plot, labels = plot_df$geneset)
  }
  
  df = res %>% 
    rename(geneset = pathway) %>%
    as_tibble() %>%
    select(-leadingEdge) %>%
    mutate(signature = input$signature_source)
  
  list(df = df, p = p)
})

# gsea results as table
output$gsea_res_table = DT::renderDataTable({
  gsea_res()$df %>%
    mutate(NES = signif(NES, 3),
           ES = signif(ES, 3),
           pval = scientific(pval),
           padj = scientific(padj)) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "none"),
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
output$mouse_hypertrophyDT = DT::renderDataTable({
  contrasts %>%
    select(modal, model, tp, MgiSymbol, logFC, FDR)%>%
    filter(!grepl("fetal", tp))%>%
    mutate(logFC = signif(logFC,3),
           FDR = scientific(FDR),
           model = as_factor(model),
           modal = as_factor(modal),
           tp = as_factor(tp)) %>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))

})

output$human_HF_bulk_indDT = DT::renderDataTable({
  contrasts_HF %>%
    mutate(AveExpr = signif(AveExpr,3),
           logFC = signif(logFC,3),
           t = signif(t,3),
           B = signif(B,3),
           P.Value = scientific(P.Value),
           adj.P.Val = scientific(adj.P.Val),
           study = as_factor(study)) %>%
    DT::datatable(escape=F, filter = "top", selection =  "multiple",
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
  contrasts %>%
    select(tp, MgiSymbol, logFC, FDR)%>%
    filter(grepl("fetal", tp))%>%
    mutate(logFC = signif(logFC,3),
           FDR = scientific(FDR),
           #model = as_factor(model),
           #modal = as_factor(modal),
           tp = as_factor(tp)) %>%
    rename(fetal_study= tp)%>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

output$human_scDT = DT::renderDataTable({
  sc.gex %>%
    mutate(Comparison= factor(Comparison),
           CellType= factor(CellType)) %>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

#### Functional analysis ####
# progeny
observeEvent(input$reset_input_TF, {
  updatePickerInput(session,
                    inputId = "select_tf", selected = character(0)
  )
})

# use switch function to update the contrast variable
contrast_ID = reactive({
  switch(input$select_contrast_func,
        # "murine_hypertrophy"= tf_hypertrophy
         "murine_hypertrophy"= list("TF" =df_tf$mm,
                                    "prog"= list("rna"= prog.res$mmRNA,
                                                 "ribo"=prog.res$mmRibo)
                                    ),
         "human_HF"= list("prog"= prog.res$hsReheat,
                          "TF" =df_tf$hs_reheat),
         "human_HF_sc"= list("prog"= prog.res$hsSC,
                             "TF" =df_tf$hs_sc),
         "human_fetal" = list("prog"= prog.res$hsfetal,
                              "TF" =df_tf$hs_fetal)
         )
         
})

# dorothea
output$tf_hypertrophy_plot = renderPlot({
    if (!is.null(input$select_tf) & ("murine_hypertrophy" %in% input$select_contrast_func))  {
      pls= map(input$select_tf, function(x){
        contrast_ID()$TF%>%
          filter(source ==x)%>%
          mutate(condition= paste(tp, modal,sep =  "_"))%>%
          ggplot(aes(x = condition, y =  score, fill = sig)) +
          facet_grid(rows= vars(model))+
          #geom_boxplot()+
          geom_hline(yintercept = 0, color= "black")+
          geom_col(width= 0.4, color ="black") +
          theme_cowplot() +
          scale_fill_manual(values = c("TRUE" = "darkgreen",
                                       "FALSE"="orange"))+
          labs(x = "experimental group", y = "TF activity score", fill = "p<0.05") +
          theme(panel.grid.major = element_line(color = "grey",
                                                size = 0.1,
                                                linetype = 1),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size= 11), 
                axis.title = element_text(size= 10)) +
          coord_flip()+
          ggtitle(x)
        
      })
      p1= cowplot::plot_grid(plotlist =  pls)
      p1
    }else if(!is.null(input$select_tf) & ("human_HF" %in% input$select_contrast_func))  {
      pls= map(input$select_tf, function(x){
        contrast_ID()$TF%>%
          #df_tf$hs_reheat%>%
          filter(source ==toupper(x))%>%
          ggplot(aes(x = study, y =  score, fill = sig)) +
          #facet_grid(rows= vars(model))+
          #geom_boxplot()+
          geom_hline(yintercept = 0, color= "black")+
          geom_col(width= 0.4, color ="black") +
          theme_cowplot() +
          scale_fill_manual(values = c("TRUE" = "darkgreen",
                                       "FALSE"="orange"))+
          labs(x = "experimental group", y = "TF activity score", fill = "p<0.05") +
          theme(panel.grid.major = element_line(color = "grey",
                                                size = 0.1,
                                                linetype = 1),
                panel.grid.minor = element_blank(),
                axis.text = element_text(size= 11), 
                axis.title = element_text(size= 10)) +
          coord_flip()+
          ggtitle(x)
        
      })
      p1= cowplot::plot_grid(plotlist =  pls)
      p1
    }else if(!is.null(input$select_tf) & ("human_HF_sc" %in% input$select_contrast_func))  {
      pls= map(input$select_tf, function(x){
        contrast_ID()$TF%>%
        #  df_tf$hs_sc%>%
          filter(source ==toupper(x))%>%
          ggplot(aes(x = celltype, y = score, fill = sig)) +
          facet_grid(rows= vars(condition))+
          geom_hline(yintercept = 0, color= "black")+
          geom_col(width= 0.4, color ="black") +
          theme_cowplot() +
          scale_fill_manual(values = c("TRUE" = "darkgreen",
                                       "FALSE"="orange"))+
          labs(x = "", y = "TF activity score", fill = "FDR<0.05") +
          theme(#panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size= 11), 
            axis.title = element_text(size= 10)) +
          coord_flip()+
          ggtitle(x)
        
      })
      p1= cowplot::plot_grid(plotlist =  pls)
      p1
    } else if(!is.null(input$select_tf) & ("human_fetal" %in% input$select_contrast_func))  {
      pls= map(input$select_tf, function(x){
        contrast_ID()$TF%>%
          #  df_tf$hs_sc%>%
          filter(source ==toupper(x))%>%
          ggplot(aes(x = condition, y = score, fill = sig)) +
          geom_hline(yintercept = 0, color= "black")+
          geom_col(width= 0.4, color ="black") +
          theme_cowplot() +
          scale_fill_manual(values = c("TRUE" = "darkgreen",
                                       "FALSE"="orange"))+
          labs(x = "", y = "TF activity score", fill = "FDR<0.05") +
          theme(#panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.text = element_text(size= 11), 
            axis.title = element_text(size= 10)) +
          coord_flip()+
          ggtitle(x)
        
      })
      p1= cowplot::plot_grid(plotlist =  pls)
      p1
    }
})

# progeny 
output$progeny_hypertropy_plot = renderPlot({
  if (input$select_contrast_func == "murine_hypertrophy") {
  
    p1 =(plot_hmap(contrast_ID()$prog$rna)+plot_hmap(contrast_ID()$prog$ribo) )
    p1
  }else if( input$select_contrast_func == "human_HF"){
    p1 =plot_hmap(contrast_ID()$prog)
    p1
  }else if( input$select_contrast_func == "human_HF_sc"){
    sc.df.prg = contrast_ID()$prog
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
  }else if( input$select_contrast_func == "human_fetal"){
    p1 =plot_hmap(contrast_ID()$prog)
    p1
  }

  })

output$progeny_table = DT::renderDataTable({
  if (input$select_contrast_func == "murine_hypertrophy") {
     lapply(names(contrast_ID()$prog), function(x){
    #df= lapply(names(test), function(x){
      contrast_ID()$prog[[x]]%>% 
      as.data.frame()%>%
      rownames_to_column("contrast")%>%
      pivot_longer(-contrast, names_to= "pathway", values_to= "score")%>% 
        mutate(sig= ifelse(abs(score)>2, T, F),
               score = signif(score, 3))
      })%>% do.call(rbind, .)%>%
      mutate(modal = factor(ifelse(grepl("rna", contrast), "rna", "ribo")),
             model = factor(ifelse(grepl("tac", contrast), "tac", "swim")),
             tp = factor(ifelse(grepl("2d", contrast), "2d", "2wk")),
             pathway= factor(pathway))%>%
      select(-contrast, )%>%
      DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
  }else  if (input$select_contrast_func == "human_HF_sc") {
    df= lapply(names(contrast_ID()$prog), function(x){
      contrast_ID()$prog[[x]] %>%
        as.data.frame()%>%
        rownames_to_column("contrast")%>%
        pivot_longer(-contrast, names_to= "pathway", values_to= "score")%>%
        mutate(celltype= factor(x),
               pathway= factor(pathway), 
               contrast= factor(contrast), 
               sig= factor(sig)
               )
    })%>% do.call(rbind, .)
    df%>%
      #prog.res$hsReheat%>%
      as.data.frame()%>%
      mutate_if(is.character, as.factor) %>% 
      mutate(sig= factor(ifelse(abs(score)>2, T, F)),
             score = signif(score, 3))%>%
      DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                    extensions = "Buttons", rownames = F,
                    option = list(scrollX = T, 
                                  autoWidth = T, 
                                  dom = "Bfrtip",
                                  buttons = c("copy", "csv", "excel")))
  }else{
  #contrast_ID()$prog%>%
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
    
  }
})

output$dorothea_table_hypertrophy = DT::renderDataTable({
  contrast_ID()$TF %>%
  #  df_tf$mm%>%
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
# 


hide("loading-content", TRUE, "fade")  

}


