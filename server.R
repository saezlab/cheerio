# SERVER
server = function(input, output, session) {
  
#### Query genes ####
# to test run some function we can create a dummy input_#input= list()
#input = list()
#input$select_gene = toupper(c("Nppb", "Nppa", "Mybpc3", "Col1a1", "Myh7", "Myh6" ))
#input$select_gene = toupper(c("GPC5", "COL4A1", "PCOLCE2", "EXT1"))
## reset the gene input button:
observeEvent(input$reset_input, {
    updatePickerInput(session,
                      inputId = "select_gene", selected = character(0)
    )
  })

## A
# cardiac mouse hypertrophy show correlation between transcript and translatome:  
output$cardiac_hyper_corr = renderPlot({
  if (!is.null(input$select_gene) ){
    
    #filter for selected genes and contrast
    to_plot_df= joint_contrast_df %>% 
      ungroup()%>%
      filter( grepl(pattern = "mm|rn", contrast_id))%>%
      select( model,modal, timepoint,gene, logFC)%>%
      filter(modal!= "proteome")%>%
      pivot_wider(names_from= modal, values_from = logFC, values_fn= mean)%>%
      mutate(labels= ifelse(gene %in% input$select_gene, gene, "background"),
             labels= factor(labels, levels= c("background", input$select_gene)),
             alphas= factor(ifelse(labels=="background", "bg","normal"))
      )%>%
      arrange(desc(labels))
    
    # assign colors!
    myColors<- myColors[1:length(levels(to_plot_df$labels))]
    names(myColors) <- levels(to_plot_df$labels)
    
    #mm models
    p1 = plot_transcipt_translat_corr(to_plot_df %>% 
                                        filter(model!= "PE"))
    
    #rat PE
    p2 = plot_transcipt_translat_corr(to_plot_df %>% 
                                        filter(model== "PE")%>%
                                        select(-timepoint))
    p2 <- remove_legend(p2)
    # add colorful facet label    
    p2= make_colorful_facet_labels(p2, "#FFE347")
    p1 = make_colorful_facet_labels(p1)
    
    # join!
    p= plot_grid(plot_grid(p2,NULL, ncol= 1),p1, ncol = 2, 
                   rel_widths = c(1,2.5),
                   labels= "AUTO")
    return(p)
  }
})

#  cardiac hypertrophy logFCs:
output$gene_expression_plots = renderPlot({
  if (!is.null(input$select_gene) )  {
    fc_vec = joint_contrast_df %>% 
      filter(gene %in% input$select_gene,
             grepl(pattern = "mm|rn", contrast_id))%>%
      pull(logFC)
    
    # separate the modal time point and model vars: 
    to_plot_df= joint_contrast_df %>% 
      filter(grepl(pattern = "mm|rn", contrast_id))%>%
      mutate(model= factor(ifelse(grepl("TAC", contrast_id), "TAC", 
                                  ifelse(grepl("swim", contrast_id ),"swim", 
                                         "PE")),levels= c("swim", "TAC", "PE")))
            
    p= map(input$select_gene, function(x){
      to_plot_df= to_plot_df %>% 
        filter(gene==x)
      
      #rows of plot df
      if(dim(to_plot_df)[1]<1){
        error_text2 <- paste(x , "\nwas not captured\nin data" )
        p <- ggplot()+ 
          annotate("text", x = 4, y = 25, size=6, label = error_text2)+
          theme_void()
        return(
          list("p"= p, "leg"= NA)
        )
      }else{
          return(
            plot_logfc_gene(red_contrast_df = to_plot_df, 
                            fg_name = "model",  
                            gene= x, 
                            max_fc = max(fc_vec), 
                            min_fc = min(fc_vec)
                            )
          )
        }
      
      
    })
    
  
    # Combine plots without legends and the legend on the right
    final_plot <- cowplot::plot_grid(
      cowplot::plot_grid(plotlist = p), # All plots without legends
      legend_lfc_plot,                 # Single legend is loaded in global.R
      ncol = 2,                                   # Arrange side by side
      rel_widths = c(length(input$select_gene)*2, 1)         # Adjust width ratio
    )
    
    return(final_plot)
    
  }
})

# Mouse heart weight: 
output$heart_weight_plot= renderPlot({
  if (!is.null(input$select_gene) )  {
    plot_hw_association(HW_DF, genes = input$select_gene)
  }
})

#IPMC_data
output$IPMC_table= DT::renderDataTable({
 ipmc_data %>%
    filter(gene %in% (input$select_gene))%>%
    select(gene, Allele, Cardiac_Hypertrophy, Other_Phenotypes,median_p_val, IPMC_link)%>%
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

## B
## human HCM bulk
output$HFgene_regulation_magnet = renderPlot({
  if (!is.null(input$select_gene)) {
    
   
    fc_vec = joint_contrast_df %>% 
      filter(contrast_id %in% hcm_contrasts,
             gene %in% input$select_gene)%>%
      pull(logFC)
    
    to_plot_df<-joint_contrast_df %>% 
      filter(contrast_id %in% hcm_contrasts)%>%
      mutate(contrast_id = factor(contrast_id))%>%
      mutate(modal= factor(modal))
    
    p= map(input$select_gene, function(x){

      to_plot_df <- to_plot_df %>%
        filter(gene== x)
        # mutate(contrast_id =factor(contrast_id, 
        #                            levels=hcm_contrasts))
      if(dim(to_plot_df)[1]<1){
        error_text2 <- paste(x , "\nwas not captured\nin data" )
        p <- ggplot()+ 
          annotate("text", x = 4, y = 25, size=6, label = error_text2)+
          theme_void()
        return(
          list("p"= p, "leg"= NA)
        )
      }else{
        
        if(length(unique(to_plot_df$modal))==2){
            fill_modal= c("#EF767A", "#8378b7")
        }else{
           fill_modal =ifelse(unique(to_plot_df$modal)=="transcriptome",
                              "#8378b7", 
                              "#EF767A")
        }
        
        return(  
          plot_logfc_gene(to_plot_df, 
                          fg_name = "modal",
                          gene= x, 
                          max_fc= max(fc_vec),
                          min_fc= min(fc_vec), 
                          fills= fill_modal
                          )
          )
      }
    })
    
    final_plot <- cowplot::plot_grid(
      cowplot::plot_grid(plotlist = p), # All plots without legends
      legend_lfc_plot,                 # Single legend is loaded in global.R
      ncol = 2,                                   # Arrange side by side
      rel_widths = c(length(input$select_gene)*2, 1)         # Adjust width ratio
    )
    return(final_plot)
  }
})

## human HCM single cell
output$HF_single = renderPlot({
  sn_contrasts <- joint_contrast_df%>%
         filter(grepl("snRNA", contrast_id))%>%
    pull(contrast_id)%>% unique()
  
  if (!is.null(input$select_gene) ) {
    to_plot_df= joint_contrast_df%>%
      filter(grepl("snRNA", contrast_id))%>%
      mutate(CellType= factor(str_extract(contrast_id, "(?<=_)[^_]+$")), 
             Comparison = factor(sapply(str_split(contrast_id, "_"), `[`, 2),
                                 levels= c( "HCMvsNF","HCMvsDCM")
             )
      )
    
    fc_vec = to_plot_df %>%
      filter(gene %in% input$select_gene)%>% 
      pull(logFC)
    
    p= map(input$select_gene, function(x){
      to_plot_df= to_plot_df%>%
        filter(gene==x)
      
      if(dim(to_plot_df)[1]<1){
        error_text2 <- paste(x , "\nwas not captured\nin data" )
        p <- ggplot()+ 
          annotate("text", x = 4, y = 25, size=6, label = error_text2)+
          theme_void()
        return(
          list("p"= p, "leg"= NA)
        )
      }else{
        return(
          plot_logfc_gene(red_contrast_df = to_plot_df, 
                          x_column =  "CellType", 
                          fg_name = "Comparison", 
                          gene= x,
                          max_fc = max(fc_vec), 
                          min_fc = min(fc_vec),
                          colored_facet= T
          )  
        )
      }  
      
    })
   
    final_plot <- cowplot::plot_grid(
      cowplot::plot_grid(plotlist = p), # All plots without legends
      legend_lfc_plot,                 # Single legend is loaded in global.R
      ncol = 2,                                   # Arrange side by side
      rel_widths = c(length(input$select_gene)*2, 1.5)         # Adjust width ratio
    )
    
    return(final_plot)
  }
  
 })

## C
## reheat (human HF): 
output$HFgene_regulation_boxplot = renderPlot({
  if (!is.null(input$select_gene) ) {
    fc_vec <- contrasts_HF %>% 
      filter(gene %in% input$select_gene)%>% 
      pull(logFC)
    
    p= map(input$select_gene, function(x){
      # prep single study df
      to_plot_df= contrasts_HF %>% 
        filter(gene== (x))%>%
        mutate( sig= adj.P.Val<0.05)%>%
        rename(contrast_id = study)
      # plot
      if(dim(to_plot_df)[1]<1){
        error_text2 <- paste(x , "was not captured\nin data" )
        return(
          ggplot() + 
            annotate("text", x = 4, y = 25, size=8, label = error_text2) + 
            theme_void()   
        )
      }else{
        return(
          plot_logfc_gene(to_plot_df,  gene= x,
                          max_fc = max(fc_vec), 
                          min_fc= min(fc_vec))
        )
      }
    })
    
    final_plot <- cowplot::plot_grid(
      cowplot::plot_grid(plotlist = p), # All plots without legends
      legend_lfc_plot,                 # Single legend is loaded in global.R
      ncol = 2,                                   # Arrange side by side
      rel_widths = c(length(input$select_gene)*2, 1)         # Adjust width ratio
    )
    
    return(final_plot)
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
                   size = 0.5, color="darkred", alpha=0.5) +
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
      geom_rug(data = sub_ranks, color= "darkred", size=1) +
      theme_classic() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      labs(x = "mean t-value", y="density")
    
    ggplotly(dens, tooltip = c("label"))
  }
})



## fetal gene expression :
output$fetal_gene_expression_plots = renderPlot({
  if (!is.null(input$select_gene) )  {
    
    fc_vec = joint_contrast_df %>% 
      filter( grepl("fetal", contrast_id ),
             gene %in% input$select_gene)%>%
      pull(logFC)
    
    p= map(input$select_gene, function(x){
      to_plot_df<- joint_contrast_df %>% 
        filter(gene==x, 
               grepl("fetal", contrast_id ))
      
      if(dim(to_plot_df)[1]<1){
        error_text2 <- paste(x , "was not captured\nin data" )
        return(
          ggplot() + 
            annotate("text", x = 4, y = 25, size=8, label = error_text2) + 
            theme_void()
          )
        }else{
          return(
            plot_logfc_gene(red_contrast_df = to_plot_df,  
                            gene= x,
                            max_fc = max(fc_vec), 
                            min_fc = min(fc_vec))
        )
      }
        
    })
    
    final_plot <- cowplot::plot_grid(
      cowplot::plot_grid(plotlist = p), # All plots without legends
      legend_lfc_plot,                 # Single legend is loaded in global.R
      ncol = 2,                                   # Arrange side by side
      rel_widths = c(length(input$select_gene)*2, 1)         # Adjust width ratio
    )
    
    return(final_plot)
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

observe({
  selected_count <- length(c(input$select_contrast_mm,
                             input$select_contrast_hs,
                             input$select_contrast_hs2
  ))
  #updateSliderInput(session, "missing_prop", max = selected_count)
  updateSliderInput(session = session, "missing_prop",
              max = selected_count,  # default, will be updated dynamically
              value = selected_count
              )
})

cont_res = eventReactive(input$submit_contrast, {
  res= get_top_consistent_gene2(joint_contrast_df = joint_contrast_df, 
                          query_contrasts = c(input$select_contrast_mm,
                                              input$select_contrast_hs,
                                              input$select_contrast_hs2),
                          cutoff = input$cut_off_genes,
                          alpha= as.numeric(input$select_alpha),
                          missing_prop = input$missing_prop
                          
                          )
  return(res)
})

output$cq_hist= renderPlot({
  cont_res()$p.hist
})

output$cq_top=renderPlot({
  cont_res()$p.top_genes
})

output$hmap_top=renderPlot({
  cont_res()$hmap_top
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

output$drugst_one_link <- renderUI({
  # Access the reactive result
  res <- cont_res()
  
  # Check if res$drugURL exists and is not NULL
  if (!is.null(res$drugst_URL)) {
    HTML(paste0(
      "<p> 3. To explore possible gene network based drug interactions, tissue expression and related disorders go to ",
      "<a href='", res$drugst_URL, "' target='_blank'> Drugst.one </a> ",
      "from the <a href='https://www.cosy.bio/' target='_blank'> Baumbachlab </a> ",
      "and explore a knowlege graph around the signature.</p>"
    ))
  } else {
    HTML("<p>No URL available. Please select contrasts and try again.</p>")
  }
})


#### Functional analysis ####
# buttons: 
observeEvent(input$reset_input_TF, {
  updatePickerInput(session,
                    inputId = "select_tf", selected = character(0)
  )
})
# input= list()
# input$select_tf = c("TRP73", "ZEB2", "AHR", "APEX1")

## PANEL Query TFS
# ## A. Animal
output$funcA_tf= renderPlot({
  if (!is.null(input$select_tf) ){
    A_contrasts<- df_func %>% 
      filter(grepl(pattern = "mm|rn", condition))%>%
      pull(condition)%>% unique()
  
    p<- get_tf_plot(contrast_ids = A_contrasts,
                  tfs= input$select_tf,
              fg_name = "model",
              colored_facet=T) 
    return(p)       
  }
 })



## B
output$funcB_tf_sc=renderPlot({
  if (!is.null(input$select_tf) ){
   sn_contrast<- df_func %>% filter(grepl("snRNA", condition))%>%
    pull(condition)%>% unique()
   
    p <- get_tf_plot(sn_contrast, 
                  tfs= input$select_tf,
                  fg_name = "model",
              colored_facet=T)
    return(p)
  }
})
output$funcB_tf_bulk=renderPlot({
  if (!is.null(input$select_tf) ){
    p<- get_tf_plot(hcm_contrasts,
              tfs= input$select_tf, 
              fg_name = "model",
              colored_facet=T)
    return(p)
  }
})

output$funcB_tf_tb_sc= DT::renderDataTable({
  df_tf$hs_sc%>%
    dplyr::select(-statistic)%>%
    mutate_if(is.character, as.factor) %>%
    mutate(score = signif(score, 3)
           #p_value = scientific(p_value)
    ) %>%
    dplyr::select(source, condition, celltype, everything())%>%
    DT::datatable(escape=F, filter = "top", selection = list(target = 'row+column'),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))
})

##C
output$funcC_tf=renderPlot({
  if (!is.null(input$select_tf) ){
   hf_contrasts<- df_func %>% 
    filter(grepl("HF|fetal|DCMvs", condition))%>%
    pull(condition)%>% unique()
  p<- get_tf_plot(hf_contrasts,
              tfs= input$select_tf,
              fg_name = NULL,
              colored_facet=FALSE)
  return(p)
  }
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


hide("loading-content", TRUE, "fade")

## PANEL QUERY CONTRASTS
observeEvent(input$reset_input_contrasts_tf, {
  updatePickerInput(session,
                    inputId = "select_contrast_hs_tf", selected = character(0)
  )
  updatePickerInput(session,
                    inputId = "select_contrast_mm_tf", selected = character(0)
  )
  updatePickerInput(session,
                    inputId = "select_contrast_hs2_tf", selected = character(0)
  )
})

observe({
  selected_count_tf <- length(c(input$select_contrast_mm_tf,
                             input$select_contrast_hs_tf,
                             input$select_contrast_hs2_tf
  ))
 
  updateSliderInput(session = session, "missing_prop_tf",
                    max = selected_count_tf,  # default, will be updated dynamically
                    value = selected_count_tf
  )
})

cont_res_tf = eventReactive(input$submit_contrast_tf, {
  res= get_consistent_tfs(df_func = df_func, 
                                query_contrasts = c(input$select_contrast_mm_tf,
                                                    input$select_contrast_hs_tf,
                                                    input$select_contrast_hs2_tf),
                                #cutoff = input$cut_off_tfs,
                                use_FDR= input$use_FDR_for_TFs, 
                                alpha= as.numeric(input$select_alpha_tf),
                                missing_prop = input$missing_prop_tf
                                
  )
  return(res)
})

output$cq_hist_tf= renderPlot({
  cont_res_tf()$p.hist
})

output$cq_top_tf=renderPlot({
  cont_res_tf()$p.top_genes
})

output$hmap_top_tf=renderPlot({
  cont_res_tf()$hmap_top
})


## PANEL Progeny
output$func_tb_pw_mm=DT::renderDataTable({
  table_to_present<- df_func%>%
    filter(database=="progeny", !grepl("hs", condition))%>%
    mutate(score= round(score, 2),
           model= factor(model),
           condition= factor(condition),
           pathway= factor(source), 
           modal = factor(ifelse(grepl("RNA", condition), "RNA", 
                                 ifelse(grepl("ribo", condition), "Ribo", "Prot")))
    )%>%
    select(-statistic, -CellType, -score_sc, -source)%>%
    select(condition, 
           modal, model, pathway , everything())
  
  
  make_nice_table(table_to_present, "score")%>%
    DT::formatSignif(columns = c("p_value", "FDR"), digits = 3)
  
  
})
output$func_tb_pw_hs=DT::renderDataTable({
  table_to_present <- df_func%>%
    filter(database=="progeny", grepl("hs", condition))%>%
    mutate(#FDR = scales::scientific(FDR, digits=3),
           #p_value= scales::scientific(p_value, 3), 
           score= round(score, 2),
           #model= factor(model),
           
           condition= factor(condition),
           pathway= factor(source), 
           comparison = factor(sapply(str_split(condition, "_"), `[`, 2)),
           modal = factor(sapply(str_split(condition, "_"), `[`, 3))
           
    )%>%
    select(-statistic, -CellType, -score_sc, -source, -model)%>%
    select(condition, comparison, 
           modal, pathway , everything())
  
  # make a nice table to plot
  make_nice_table(table_to_present, "score")%>%
    DT::formatSignif(columns = c("p_value", "FDR"), digits = 3)
  
  
})



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
# sig = joint_contrast_df %>% filter(cc == "A")
# gss = example_geneset

# perform GSEA with plotting
gsea_res = eventReactive(input$submit, {
  if (ncol(gs()) == 1) {
    df= map(unique(signature()$contrast_id), function(x){
      y= signature()%>% filter(contrast_id == x)
      vect= deframe(y[,c("gene", "logFC")])
      res = fgsea(stats= vect ,
                  pathways= gs()$gene, 
                  nperm = 1000)
      as.data.frame(res)%>% 
        mutate(cc= x)
    }) %>% do.call(rbind,. )
    
    p1= df %>% 
      ggplot(., aes(x= cc, y= NES, fill = padj<0.05)) +
      geom_col(width= 0.4, color ="black") +
      scale_fill_manual(values = c("TRUE" = "#4D7298",
                                   "FALSE"="grey", 
                                   drop= FALSE))+
      coord_flip()+
      labs(x= "contrast ID", 
           y= "normalized enrichment score (NES)", 
           fill = "FDR<0.05")+
      theme(panel.grid.major = element_line(color = "grey",
                                            linewidth = 0.1,
                                            linetype = 1),
            panel.border = element_rect(fill= NA, linewidth=1, color= "black"), 
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.text = element_text(size= 11), 
            axis.title = element_text(size= 10))
    
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
      facet_grid(rows= vars(pathway))+
      geom_col(width= 0.4, color ="black") +
      scale_fill_manual(values = c("TRUE" = "#4D7298",
                                   "FALSE"="grey", 
                                   drop= FALSE))+
      coord_flip()+
      labs(x= "contrast ID", 
           y= "normalized enrichment score (NES)", 
           fill = "FDR<0.05")+
      theme(panel.grid.major = element_line(color = "grey",
                                            linewidth = 0.1,
                                            linetype = 1),
            panel.border = element_rect(fill= NA, linewidth=1, color= "black"), 
            panel.grid.minor = element_blank(),
            axis.line = element_blank(),
            axis.text = element_text(size= 11), 
            axis.title = element_text(size= 10))
    
   
  }
  
df = df %>% 
  dplyr::select(cc, everything())%>%
  dplyr::rename(geneset = pathway,
                contrast_id = cc) %>%
    as_tibble() %>%
    #select(-leadingEdge) %>%
    mutate(signature = input$signature_source)
  
  list(df = df, p = p1)
})

# gsea results as table
output$gsea_res_table = DT::renderDataTable({
  gsea_res()$df%>%
    mutate(NES = signif(NES, 3),
           ES = signif(ES, 3),
           geneset= factor(geneset), 
           contrast_id = factor(contrast_id))%>%
           # pval = scientific(pval),
           # padj = scientific(padj)) %>%
    DT::datatable(escape=F, filter = "top", selection = list(target = "column"),
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))%>%
    DT::formatSignif(columns = c("pval", "padj"), digits = 3)
})

# gsea results as plots
output$gsea_res_plots = renderPlot({
  gsea_res()$p
})



####Download center ####

# make gene contrast data accessible:
output$mouse_hypertrophyDT = DT::renderDataTable( {
  
  joint_contrast_df%>%
    filter(cc== "A")%>%
    select( model,modal, timepoint,gene, logFC, FDR_mod)%>%
    mutate(logFC = signif(logFC,1),
           FDR = FDR_mod
           ) %>%
      select(-FDR_mod)%>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T, 
                                autoWidth = T, 
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))%>%
    DT::formatSignif(columns = c("FDR"), digits = 3)
  
})


output$human_HF_bulk_summDT = DT::renderDataTable({
  
  ranks %>%
    mutate(mean_lfc = signif(mean_lfc,3),
           mean_t = signif(mean_t,3)) %>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))%>%
  DT::formatSignif(columns = c("fisher_pvalue"), digits = 3)
  })

output$human_fetalDT = DT::renderDataTable({
  joint_contrast_df%>%
    select(contrast_id, gene, logFC, FDR,sig)%>%
    filter(grepl("fetal", contrast_id))%>%
    mutate(logFC = signif(logFC,3),
           #FDR = scientific(FDR),
           contrast_id = factor(contrast_id)) %>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                  extensions = "Buttons", rownames = F,
                  option = list(scrollX = T,
                                autoWidth = T,
                                dom = "Bfrtip",
                                buttons = c("copy", "csv", "excel")))%>%
    DT::formatSignif(columns = c("FDR"), digits = 3)
  
  
})

output$human_HCMDT = DT::renderDataTable({
  
  df<- joint_contrast_df%>% 
    filter(grepl(pattern="HCM", contrast_id) )%>%
    #mutate(modality= factor(ifelse(grepl("bulk", contrast_id),"bulk", "single_cell")))%>%
    mutate(logFC = signif(logFC,3),
         #FDR = scientific(FDR),
         contrast_id = factor(contrast_id)) %>%
    mutate(CellType= factor(str_extract(contrast_id, "(?<=_)[^_]+$")), 
           Comparison = factor(sapply(str_split(contrast_id, "_"), `[`, 2),
                               levels= c( "HCMvsNF","HCMvsDCM")
           ))%>%
    select(contrast_id, modal, CellType,Comparison, gene, logFC, FDR,sig)%>%
    DT::datatable(escape=F, filter = "top", selection = "multiple",
                                 extensions = "Buttons", rownames = F,
                                 option = list(scrollX = T,
                                               autoWidth = T,
                                               dom = "Bfrtip",
                                               buttons = c("copy", "csv", "excel")))%>%
    DT::formatSignif(columns = c("FDR"), digits = 3)
    

})


}


