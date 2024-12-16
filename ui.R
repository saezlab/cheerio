# UI
source("sub/global.R")
source("sub/helper.R")
ui = function(request) {
  fluidPage(
    rclipboardSetup(),
    #tags$script(src = "https://kit.fontawesome.com/acea36f561.js"),
    useShinyjs(),
    # HTML("<script src='https://cdn.drugst.one/latest/drugstone.js'></script>",
    #      "<link rel='stylesheet' href='https://cdn.drugst.one/latest/styles.css'>"),
    
    navbarPage(
      id = "menu", 
      theme = shinytheme("flatly"),
     # theme = bs_theme(preset = "Simplex"),
      title = div(img(src="logo_cheerio-removebg-preview.png", width="25", height="25"),
                  "CHEERIO"),
      windowTitle = "CHEERIO",
      collapsible=T,
      
      #### Welcome ####
      tabPanel(
        title = "Welcome",
        icon = icon("home"),
        sidebarPanel(
          includeMarkdown("inst/landingpage_sidebar.md")
        ),
        mainPanel(
          includeMarkdown("inst/landingpage.md"),
          DT::dataTableOutput("meta_table")%>%
            withSpinner()
        )
      ),
      
      #### Query Genes ####
      tabPanel(
        title = "Query Genes",
        icon = icon("search"),
        sidebarPanel(
          includeMarkdown("inst/query_genes_sidebar.md"),
          pickerInput(inputId = "select_gene", 
                      label = "Select gene(s)",
                      choices = (sort(unique(joint_contrast_df$gene))), 
                      multiple = T,
                      options = list(`live-search` = TRUE,
                                     size=10, `max-options` = 6)
                     # selected = toupper(c("Nppb", "Nppa", "Mybpc3", "Col1a1", "Myh7", "Myh6" )) 
                      ),
          
          actionButton("reset_input", "Reset genes")
  
        ),
        mainPanel(
          tabsetPanel(
            type = "tabs",
            tabPanel("A. Animal models",
                     h3("1. Regulation in rodent cardiac hypertrophy models"),
                     h4("1.1 Gene expression regulation"),
                     plotOutput("gene_expression_plots", width = "100%", height = "600px")%>%
                       withSpinner(),#, width = "100%", height = "600px"),
                     p("Bar plots indicate log2 fold changes for different contrasts and data sets. Please refer to the table at the landing page for 
                     more information on the contrasts."),
                     p(strong("Abbreviations:")),
                     p("mm: mouse; rn: rat; RNA: RNA-seq (transcriptional regulation); ribo: Ribo-seq (translational regulation); 2d: two days; 2wk: two weeks; 
                     swim: swimming (physiologic hypertrophy); TAC: transverse-aortic-constriction (pathologic hypertrophy); PE: Phenylephrine"),
                     
                     br(),
                     br(),
                     h4("1.2  Ribo-seq and RNA-seq correlation"),
                     plotOutput("cardiac_hyper_corr")%>%
                       withSpinner(),
                     p("Comparison of logFC values in Ribo-seq (y-axis) and RNA-seq (x-axis) for different models (A - in vitro, B - in vivo)"),
                     p(strong("Abbreviations:")),
                     p("RNA: RNA-seq (transcriptional regulation); ribo: Ribo-seq (translational regulation); 2d: two days; 2wk: two weeks; 
                     swim: swimming (physiologic hypertrophy); TAC: transverse-aortic-constriction (pathologic hypertrophy); PE: Phenylephrine"),
                     
                     br(),
                     br(),
                     br(),
                     hr(),
                     
                     #hw
                     h3("2. Phenotype associations"),
                     h4("2.1 Gene expression association with heart weight in mice after swim or TAC"),
                     br(),
                     p(br()),
                     plotOutput("heart_weight_plot", width = "100%", height = "400px")%>%
                       withSpinner(),
                     p("Univariate Linear models are used to describe the association between gene expression and normalized heart weight.
                        The model (normalized heart weight ~ b0 + b1 gene expression) was fit for each gene (panels) and contrast (x-axis).
                       The coefficient of gene expression is displayed (y-axis), while the shape indicates whether the p-value of the coefficent was significant (p-value <0.05, triangle) or not (p-value > 0.05, circle). 
                       The color indicates the R² of the model."), 
                     br(),
                     ##ipmc table
                    h4("2.2 Known associations of genetic variants with mouse phenotypes"),
                    h5("The IMPC (International Mouse Phenotyping Consortium) gathers data on genetic mouse models and a diversity of measured phenotypes. We show whether the selected 
                       genes are associated with any cardiovascular phenotypes in the IMPC database. For more information, click on the gene name for a direct link to the corresponding IMPC entry or on the Jax link for the MGI site"),
                    br(),
                    DT::dataTableOutput("IPMC_table")%>%
                      withSpinner(),
                   hr()
                     ),
            tabPanel("B. Human Cardiac Hypertrophy", 
                     ##Magnet
                     h4("Regulation on bulk level"),
                     
                     plotOutput("HFgene_regulation_magnet", width = "100%", height = "500px")%>%
                       withSpinner(),
                     p("hs: human; HCM: hypertrophic cardiomyopathy; HCMrEF: hypertrophic cardiomyopathy with reduced ejection fraction; HCMpEF: hypertrophic cardiomyopathy with preserved ejection fraction; cHYP: compensated cardiac hypertrophy (non-failing); DCM: dilated cardiomyopathy; NF: non-failing healthy heart"),
                     br(),
                     br(),
                     hr(),
                     
                     h4("Regulation on single cell level (human)"),
                     #h6("Celltype regulation compared between non failing, DCM and HCM"),
                     
                     plotOutput("HF_single", height= "900px", width= "100%") %>%
                       withSpinner(),
                     p("HCM: hypertrophic cardiomyopathy; DCM: dilated cardiomyopathy; NF: non-failing healthy heart; VSMC: Vascular smooth muscle cell; PC: Pericyte; Neu: Neuronal; MP: Macrophage; MC: Mast cell; Lympho: Lymphocyte; LEC: Lymphatic endothelial; FB: Fibroblast; Endo: Endocard ; ECl: Endothelial; CM: Cardiomyocyte; Adipo: Adipocyte"),
                     
                     br(),
                     br(),
                     br(),
                     hr()
            ),
            tabPanel("C. Human Heart Failure",
                     h4("Regulation in bulk transcriptomics from human heart failure studies (DCM and ICM patients)"),
                     h6("Heart failure bulk transcriptomic studies"),
                     plotOutput("HFgene_regulation_boxplot", width = "100%", height = "500px") %>%
                       withSpinner(),
                     
                     h5("Distribution of mean t-values"),
                     h6("Mean t-values distribution of all bulk transcriptomics studies."),
                     plotlyOutput("mean_t_dist", width = "100%", height = "250px") %>%
                       withSpinner(),
                     p("Position of a gene towards the  left indicates strong and significant downregulation across all studies compared to all other genes."),
                     p("Position of a gene towards the right indicates strong and significant upregulation across all studies compared to all other genes."),
                     
                     h5("Ranking of queried genes"),
                     h6("Consensus ranking, the lower the rank the more consistently is the gene significantly regulated in human HF"),
                     plotlyOutput("rank_position", width = "100%", height = "125px") %>%
                       withSpinner(),
                     
                     br(),
                     br(),
                     br(),
                     hr(),
                     
                    
                     ),
            tabPanel("D. Fetal gene program", 
                     h3("Regulation in Fetal vs. Adult Human Hearts"),
                     plotOutput("fetal_gene_expression_plots", height= "250px"),
                     p(strong("Abbreviations:")),
                     p("hs: human; RNA: RNA-Seq"),
                     p("a positive logFC indicates that the gene is higher expressed in the fetal human heart compared to the adult human heart"),
                     p("a negative logFC indicates that the gene is lower expressed in the fetal human heart compared to the adult human heart"),
                     
                     hr()
                     )
            )
          )
    
      ),
      
      

# Query contrasts -----------------------------------------------------------------------------
tabPanel(
  title = "Custom Signature",
  icon = icon("compress-alt"),
  sidebarPanel(
    includeMarkdown("inst/query_contrasts_sidebar.md"),
    br(),
    #hr(),
    strong("1. Select contrast(s)"),
    br(), 
    pickerInput(inputId = "select_contrast_mm", 
                label = "A) Animal studies:",
                choices = joint_contrast_df%>% 
                  filter(cc =="A")%>% 
                  arrange(contrast_id)%>%
                  pull(contrast_id)%>%
                  unique(), 
                multiple = T,
                options = list(`live-search` = TRUE,
                               size=10, `max-options` = 10),
                # selected = c("mm_TAC_RNA_2w") 
    ),
    
    pickerInput(inputId = "select_contrast_hs", 
                label = "B) Human HCM studies:",
                choices = joint_contrast_df%>% 
                  filter(cc =="B")%>%
                  arrange(contrast_id)%>%
                  pull(contrast_id)%>% 
                  unique(), 
                multiple = T,
                options = list(`live-search` = TRUE,
                               size=10, `max-options` = 10),
                # selected = c("hs_HCMvsNF_snRNA_CM", "hs_HCMvsNF_RNA") 
    ),
    
    pickerInput(inputId = "select_contrast_hs2", 
                label = "C+D) Human HF or fetal gene program",
                choices = joint_contrast_df%>% 
                  filter(cc %in% c("C", "D"))%>%
                  arrange(contrast_id)%>%
                  pull(contrast_id)%>% unique(), 
                multiple = T,
                options = list(`live-search` = TRUE,
                               size=10, `max-options` = 10),
                # selected = c("hs_fetal_RNA") 
    ),
   
    actionButton("reset_input_contrasts", "Reset contrasts"),
    br(),
    #br(),
    hr(),
    pickerInput(inputId = "select_alpha", 
                label = "2. Select alpha for FDR cut off:",
                choices = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2), 
                multiple = F,
                options = list(`live-search` = TRUE,
                               size=10, `max-options` = 10),
                selected = 0.05
                ),
    
    sliderInput("cut_off_genes", "3. Select number of genes to plot:",
                min = 1, max = 70,
                value = 10, step= 1),
    hr(),
    sliderInput("missing_prop", "4. Select minimum number of contrasts to report a gene",
                min = 1, max = 10,
                value = 1, step= 1),
    br(), 
    strong("5. Compare contrasts!"), 
    br(), 
    actionButton("submit_contrast", label="Submit",
                 icon=icon("paper-plane")) 
  ),
  mainPanel(
    h3("1. Search for consistent genes"),
    h4("1.1 Dysregulated genes (FDR-cutoff)"),
    plotOutput("cq_hist")%>%withSpinner(),
    p("A. Histogram showing the number of shared differntially expressed genes (DEGs) (y-axis) reported by differnt number of contrasts (x-axis).
                       Blue colored bars indicate which DEGs will be considered based on the selected minimum number of contrasts where a gene is significantly regulated."),
    p("B. The consistency of direction of regulation for the selected DEGs. Only consistent DEGs will be further considered. A gene is considered consistent if all datasets where a gene is significantly regulated show the same directionality (up- or downregulated)."),
    hr(),
    br(),
    br(),
    h4("1.2 Full signature of consistently regulated genes"),
    
    plotOutput("hmap_top")%>% withSpinner(),
    p("Heatmap displaying the full signature of consistently regulated genes. A black square indicates that the gene was not detected in the respective dataset."),
    hr(),
    br(),
    br(),
    
    h4("1.3 Top consistently up- and downregulated genes"),
    plotOutput("cq_top")%>%withSpinner(),
    p("Boxplots display top upregulated (top panel) and downregulated (bottom panel) genes from the generated signature."),
    hr(),
    br(),
    br(),
    h3("2. Characterize the shared gene signature"),
    br(),
    h5("You can now explore functional annotations of the up or downregulated genes!"),
    p("1.1 Click one of the buttons to copy the up- or downregulated genes into your clipboard."),
       br(),
    fluidRow(
      uiOutput("clipup",  style = 'display: inline-block;margin-left: 15px;' ),
    #  br(),
      uiOutput("clipdn",  style = 'display: inline-block' )
      ),
    br(),
    # HTML("<p> 2a. To check the selected genes for footprints of TFs or Pathways go to <a href='https://saezlab-funki-analysis-rnfy3j.streamlit.app/'> FUNKi </a> 
    #     from the <a href=' https://saezlab.org'> Saezlab </a and paste the genes into the gene submission field. 
    #     </p> "),
    HTML("<p> 1.2 To enrich gene sets from various biological databases (GO, MSIG, KEGG etc.) go to <a href='https://maayanlab.cloud/Enrichr/'> Enrichr </a> 
        or <a href='https://maayanlab.cloud/enrichr-kg'> Enrichr-KG </a> from the <a href=' https://labs.icahn.mssm.edu/maayanlab/'> Mayan lab</a and paste the genes into the gene submission field. 
        </p> "),
    uiOutput("drugst_one_link"),
    p("3. Explore possible functional processes that your selected contrasts have in common and generate a hypothesis!"),
    hr(),  
    br(),
    br(),
    h3("3. Explore the signature in table form"),
    tabsetPanel(
      type = "tabs",
      tabPanel("Upregulated",
               DT::dataTableOutput("cq_table_up")),
      tabPanel("Downregulated",
               DT::dataTableOutput("cq_table_dn")),
      tabPanel("all", DT::dataTableOutput("cq_table_full"))
    )%>%withSpinner(),
    hr(),  
    br(),
    br(),
    
     
  )
),

      
      #### Functional analysis ####
tabPanel(
        title = "Functional Genomics",
        icon = icon("chart-line"),
        tabsetPanel(
          type = "tabs",
          tabPanel("Query a TF", 
                   br(),
            sidebarPanel(
              includeMarkdown("inst/functional_analysis_sidebar.md"),
              p(""),
              pickerInput(inputId = "select_tf",
                      label = "Select TF(s)",
                      choices = sort(TFs),
                      multiple = T,
                      options = list(`live-search` = TRUE,
                                     size=6, `max-options` = 6),
                      # selected = c("GATA4", "GATA3", "HIF1A", "HIF1B")
                      ),
              actionButton("reset_input_TF", "Reset TFs")
              ),
            mainPanel(
              tabsetPanel(
                type = "tabs",
                tabPanel("A. Animal models",
                         h3("Transcription factor activities"),
                         plotOutput("funcA_tf", width = "100%", height= "800px"),#, width = "100%", height = "600px"),
                         p("Bar plots indicate TF activties (ulm scores, see decoupleR) for different contrasts and 
                         data sets. Please refer to the table at the landing page for more information on the contrasts."),
                         p(strong("Abbreviations:")),
                         p("mm: mouse; rn: rat; RNA: RNA-seq (transcriptional regulation); ribo: Ribo-seq (translational regulation); 2d: two days; 2wk: two weeks; 
                     swim: swimming (physiologic hypertrophy); TAC: transverse-aortic-constriction (pathologic hypertrophy); PE: Phenylephrine"),
        
                         br(),
                         br(),
                         hr()
                         ),
                tabPanel("B. Human Hypertrophic Cardiomyopathy",
                         ##Magnet
                         h3("Estimated transcription factor activities - bulk transcriptomics HCM"),
                         plotOutput("funcB_tf_bulk", width = "100%", height = "500px"),
                         p("Bar plots indicate TF activties (ulm scores, see decoupleR) for different contrasts and 
                         data sets. Please refer to the table at the landing page for more information on the contrasts."),
                         p(strong("Abbreviations:")),
                         p("hs: human; HCM: hypertrophic cardiomyopathy; HCMrEF: hypertrophic cardiomyopathy with reduced ejection fraction; HCMpEF: hypertrophic cardiomyopathy with preserved ejection fraction; cHYP: compensated cardiac hypertrophy (non-failing); DCM: dilated cardiomyopathy; NF: non-failing healthy heart"),
                         br(),
                         br(),
                         hr(),
                         h3("Estimated transcription factor activities - single cell transcriptomics HCM"),
                         plotOutput("funcB_tf_sc", width = "100%", height = "500px"),
                         p("Bar plots indicate TF activties (ulm scores, see decoupleR) for different contrasts and 
                         data sets. Please refer to the table at the landing page for more information on the contrasts."),
                         p(strong("Abbreviations:")),
                         p("HCM: hypertrophic cardiomyopathy; DCM: dilated cardiomyopathy; NF: non-failing healthy heart; VSMC: Vascular smooth muscle cell; PC: Pericyte; Neu: Neuronal; MP: Macrophage; MC: Mast cell; Lympho: Lymphocyte; LEC: Lymphatic endothelial; FB: Fibroblast; Endo: Endocard ; ECl: Endothelial; CM: Cardiomyocyte; Adipo: Adipocyte"),
                         br(),
                         br(),
                         hr()
                         ),
                tabPanel("C+ D. Human Heart Failure and Fetal Program",
                         h3("Transcription factor activities"),
                        plotOutput("funcC_tf", width = "100%", height = "500px") %>%
                         withSpinner(),
                         p("Bar plots indicate TF activties (ulm scores, see decoupleR) for different contrasts and 
                         data sets. Please refer to the table at the landing page for more information on the contrasts."),
                         p("Fetal: a positive logFC indicates that the gene is higher expressed in the fetal human heart compared to the adult human heart, a negative logFC indicates that the gene is lower expressed in the fetal human heart compared to the adult human heart"),
                         p(strong("Abbreviations:")),
                         p("hs: human; HF: heart failure; NF: non-failing healthy heart; RNA: RNA-Seq"),
                        br(),
                        br(),
                        hr()
                        )
                )# tabset panel (ABCD)
              )# main panel ()
            ), # query TF panel
          tabPanel("Conserved TFs", 
                   br(),
                   sidebarPanel(
                     includeMarkdown("inst/conservedTFs_sidebar.md"),
                     br(),
                     #hr(),
                     strong("1. Select contrast(s)"),
                     br(), 
                     pickerInput(inputId = "select_contrast_mm_tf", 
                                 label = "A) Animal studies:",
                                 choices = sort(df_func%>% filter(cc =="A")%>% pull(condition)%>% unique()), 
                                 multiple = T,
                                 options = list(`live-search` = TRUE,
                                                size=10, `max-options` = 10),
                                #selected = c("mm_TAC_RNA_2w", "mm_TAC_ribo_2w") 
                     ),
                     
                     pickerInput(inputId = "select_contrast_hs_tf", 
                                 label = "B) Human HCM studies:",
                                 choices = sort(df_func%>% filter(cc =="B")%>% pull(condition)%>% unique()), 
                                 multiple = T,
                                 options = list(`live-search` = TRUE,
                                                size=10, `max-options` = 10)
                                 #selected = c("hs_HCMvsNF_snRNA_CM", "hs_HCMvsNF_RNA") 
                     ),
                     
                     pickerInput(inputId = "select_contrast_hs2_tf", 
                                 label = "C+D) Human HF or fetal gene program",
                                 choices = df_func%>% filter(cc %in% c("C", "D"))%>% pull(condition)%>% unique(), 
                                 multiple = T,
                                 options = list(`live-search` = TRUE,
                                                size=10, `max-options` = 10)
                                 #selected = c("hs_fetal_RNA") 
                     ),
                     
                     actionButton(inputId = "reset_input_contrasts_tf", label = "Reset contrasts"),
                     br(),
                     #br(),
                     hr(),
                     pickerInput(inputId = "select_alpha_tf", 
                                 label = "2.1 Select alpha level:",
                                 choices = c(0.0001, 0.001, 0.01, 0.05, 0.1, 0.2), 
                                 multiple = F,
                                 options = list(`live-search` = TRUE,
                                                size=10, `max-options` = 10),
                                 selected = 0.05
                     ),
                     br(),
                     strong("2.2 Use multiple hypothesis correction?"), 
                    
                     switchInput(inputId = "use_FDR_for_TFs",
                                 onLabel = "Yes", offLabel = "No", value=T, 
                                 size= "mini"),
                     p("Caution! Not using FDR correction will produce more false positive results and should only be used for hypothesis generation!"), 
                     # sliderInput("cut_off_tfs", "3. Select number of TFs to plot:",
                     #             min = 1, max = 70,
                     #             value = 10, step= 1),
                     hr(),
                     sliderInput("missing_prop_tf", "3. Select minimum number of contrasts to report a TF",
                                 min = 1, max = 10,
                                 value = 1, step= 1),
                     br(), 
                     strong("4. Get conserved TFs!"), 
                     br(), 
                     actionButton(inputId = "submit_contrast_tf", label="Submit",
                                  icon=icon("paper-plane")) 
                   ),
                   mainPanel(
                     h3("Search for conserved TFs"),
                     h4("Active TFs (FDR-cutoff)"),
                     plotOutput("cq_hist_tf")%>%withSpinner(),
                    p("A. Histogram showing the number of shared estimated transcription factors activity (y-axis) reported by differnt number of contrasts (x-axis).
                       Blue colored bars indicate which transcription factors will be considered based on the selected minimum number of contrasts."),
                     p("B. The consistency in direction of activity is shown for the considered TFs. Only consistent TFs will be further considered. A TF is considered consistent if all datasets where a TF is estimated to be significantly regulated show the same directionality (up- or downregulated)."),
                     p(""),
                     hr(),
                     br(),
                     br(),
                     h4("Full panel of TF(s)"),
                     plotOutput("hmap_top_tf")%>%withSpinner(),
                     p("Heatmap showing estimated TF activities across selected contrasts. Black indicates missing values."),
                     hr(),
                     br(),
                     br(),
                     
                     h4("Top conserved TFs"),
                     plotOutput("cq_top_tf")%>%withSpinner(),
                     p("Boxplot with jitter points displays top conserved TFs based on mean activity scores."),
                     hr(),
                     br(),
                     br()
                     )#main panel conserved tfs
                   ),# conserved TF panel
          tabPanel("Pathway activities", 
                   br(),
                   sidebarPanel(
                     includeMarkdown("inst/pathways_sidebar.md")
                   ),
                   mainPanel(                    
                     h3("Pathway activities"),
                     tabsetPanel(
                       type = "tabs",
                       tabPanel("Animal Models",
                                DT::dataTableOutput("func_tb_pw_mm")),
                       tabPanel("Humans",
                                DT::dataTableOutput("func_tb_pw_hs"))                            
                       )
                     )
                   
          ),
          )# TF- pway- tf search
        ),

      #### Input data ####
      tabPanel(
        title = "Enrichment Analysis",
        icon = icon("file-upload"),
        sidebarPanel(
          includeMarkdown("inst/input_data_sidebar.md"),
          h5("Use example data?"),
          switchInput(inputId = "take_example_data",
                      onLabel = "Yes", offLabel = "No", value=F), 
          p("Gene sets must be uploaded as .csv file. The gene set members must 
             be stored in column named 'gene'. In case of multiple gene sets a 
             second column named 'geneset' must be added containing the gene set 
             name/identifier."),
          fileInput("user_input", label="Upload gene sets (.csv)"),
          pickerInput(inputId = "select_contrast_enrichment", 
                      label = "Contrast(s)",
                      choices = c("A. Animal models", "B. Human HCM", "C. Human HF", "D. Fetal reprogramming"), 
                      multiple = F,
                      selected = c("A. Animal models") 
          ),
          p("GSEA will be performed upon clicking the submit button."),
          actionButton("submit", label="Submit",
                       icon=icon("paper-plane")) 
          ),
        mainPanel(
          h4("GSEA plots"),
          plotOutput("gsea_res_plots", width = "100%", height = "600px") %>%
            withSpinner(),
          hr(),
          h4("GSEA result"),
          DT::dataTableOutput("gsea_res_table") %>%
            withSpinner()
        )
      ),
      
      #### Download Center ####
      tabPanel(
        title = "Access Data",
        icon = icon("table"),
        sidebarPanel(
          includeMarkdown("inst/meta_analysis_results_sidebar.md")
        ),
        mainPanel(
          h4("Query data in a convenient table format"),
          h6("Select your desired contrasts, sort and filter data in table format"),
          tabsetPanel(
            type = "tabs",
            tabPanel("A. Animal Models", DT::dataTableOutput("mouse_hypertrophyDT"),
                     br(), p("Download full data from zenodo, LINK")),
            tabPanel("B. Human HCM", DT::dataTableOutput("human_HCMDT"),
                     br(), p("Download full data from zenodo, LINK")),
           tabPanel("C. Human HF",  DT::dataTableOutput("human_HF_bulk_summDT"),
                     br(), p("Download full data from zenodo, LINK")),
            tabPanel("D. Fetal gene program", DT::dataTableOutput("human_fetalDT"),
                     br()
                     )
           ),
          br(), 
          hr(),
          h4("How to access full data"),
          br(),
          p("We provide a full data set of contrasts to be downloaded
            from Zenodo. LINK"),
          p("You can download the .zip file and access the processed data")
          
          
        )
      )
      
      #### Footer ####
      #footer = column(12, align="center", "CHEERIO-App")
    ) # close navbarPage
  ) # close fluidPage
}
