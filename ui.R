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
      title = div(img(src="logo_saezlab.png", width="25", height="25"),
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
          includeMarkdown("inst/landingpage.md")
        )
      ),
      
      #### Query genes ####
      tabPanel(
        title = "Query genes",
        icon = icon("search"),
        sidebarPanel(
          includeMarkdown("inst/query_genes_sidebar.md"),
          pickerInput(inputId = "select_gene", 
                      label = "Select gene(s)",
                      choices = (sort(unique(joint_contrast_df$gene))), 
                      multiple = T,
                      options = list(`live-search` = TRUE,
                                     size=10, `max-options` = 6),
                      selected = toupper(c("Nppb", "Nppa", "Mybpc3", "Col1a1", "Myh7", "Myh6" )) 
                      ),
          
          actionButton("reset_input", "Reset genes")
        #   
        #   checkboxGroupInput(inputId= "contrasts", label= "Select contrasts", 
        #                      selected = c("Murine_Hypertrophy" , "Human_HF", "Fetal" ),
        #                      inline = FALSE, width = NULL, 
        #                      choiceNames = c("Murine Hypertrophy" , "Human heart failure", "Human fetal" ),
        #                      choiceValues = c("Murine_Hypertrophy" , "Human_HF", "Fetal" ))
        ),
        mainPanel(
          tabsetPanel(
            type = "tabs",
            tabPanel("A. Animal models",
                     h3("Regulation in Murine Cardiac Hypertrophy models"),
                     h4("Gene expression regulation"),
                     plotOutput("gene_expression_plots", width = "100%", height = "600px")%>%
                       withSpinner(),#, width = "100%", height = "600px"),
                     p(strong("Abbreviations:")),
                     p("ribo, Ribo-seq (translational regulation); rna, RNA-seq (transcriptional regulation); 2d, two days; 2wk, two weeks; 
               swim, swimming (physiologic hypertrophy); tac, transverse-aortic-constriction (pahtologic hypertrophy)"),
                     
                     br(),
                     br(),
                     h4("Ribo-seq and RNA-seq correlation"),
                     plotOutput("cardiac_hyper_corr")%>%
                       withSpinner(),
                     p(strong("Abbreviations:")),
                     p("2d, two days; 2wk, two weeks; swim, swimming (physiologic hypertrophy); tac, transverse-aortic-constriction (pahtologic hypertrophy)"),
                     
                     br(),
                     br(),
                     br(),
                     hr(),
                     
                     #hw
                     h3("Phenotype associations"),
                     h4("Explore posible associations of a gene's expression with heart weight"),
                     br(),
                     p(br()),
                     plotOutput("heart_weight_plot", width = "100%", height = "400px")%>%
                       withSpinner(),
                     br(),
                     ##ipmc table
                    h4("Explore posible associations of genetic variants with mouse phenotypes"),
                    h5("Genes are queried in IMPC (International Mouse Phenotyping Consortium"),
                    br(),
                    DT::dataTableOutput("IPMC_table")%>%
                      withSpinner(),
                   hr()
                     ),
            tabPanel("B. Human Hypertrophic Cardiomyopathy", 
                     ##Magnet
                     h4("Regulation on bulk level"),
                     
                     plotOutput("HFgene_regulation_magnet", width = "100%", height = "500px")%>%
                       withSpinner(),
                     p("HCM, hypertrophic cardiomyopathy; DCM, dilated cardiomyopathy; NF, non-failing"),
                     br(),
                     br(),
                     hr(),
                     
                     h4("Regulation on single cell level"),
                     #h6("Celltype regulation compared between non failing, DCM and HCM"),
                     
                     plotOutput("HF_single", height= "900px", width= "100%") %>%
                       withSpinner(),
                     p("HCM, hypertrophic cardiomyopathy; DCM, dilated cardiomyopathy; NF, non-failing"),
                     
                     br(),
                     br(),
                     br(),
                     hr()
            ),
            tabPanel("C. Human Heart Failure",
                     h4("Regulation in bulk from human heart failure studies"),
                     h6("HF bulk transcripomic studies"),
                     plotOutput("HFgene_regulation_boxplot", width = "100%", height = "500px") %>%
                       withSpinner(),
                     
                     h5("Distribution of mean t-values"),
                     h6("Mean t-values distribution of all bulk studies"),
                     plotlyOutput("mean_t_dist", width = "100%", height = "250px") %>%
                       withSpinner(),
                     
                     h5("Ranking of queried genes"),
                     h6("Consensus ranking, the lower the rank the more consistently is the gene regulated in human HF"),
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
                     p("rna_fetal1, fetal study 1 (Spurell et al, 2022); rna_fetal2, fetal study 2 (Akat et al, 2014)"),
                     
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
                selected = c("mm_TAC_RNA_2w") 
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
                selected = c("hs_HCMvsNF_snRNA_CM", "hs_HCMvsNF_RNA") 
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
                selected = c("hs_fetal_RNA") 
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
    h3("Search for consistent genes"),
    h4("Dysregulated genes (FDR-cutoff)"),
    plotOutput("cq_hist")%>%withSpinner(),
    p("A. Barplot displays the overlap between the selected contrasts. Red indicates the number of genes for Signature generation."),
    p("B. Barplot displays the number of genes of the intersect left that and their directionality of the intersection of all selected contrasts, i.e. how many genes are commonly up- or down- or inconsistently (up and down) regulated."),
    p(""),
    hr(),
    br(),
    br(),
    h4("Full signature of dysregulated gene(s)"),
    
    plotOutput("hmap_top")%>%withSpinner(),
    hr(),
    br(),
    br(),
    
    h4("Top consistently dysregulated gene(s)"),
    plotOutput("cq_top")%>%withSpinner(),
    hr(),
    br(),
    br(),
    h3("Characterize the shared gene signature"),
    br(),
    h5("You can now explore functional annotations of the up or downregulated genes!"),
    p("1. Click one of the buttons to copy genes into your clipboard."),
       br(),
    fluidRow(
      uiOutput("clipup",  style = 'display: inline-block;margin-left: 15px;' ),
    #  br(),
      uiOutput("clipdn",  style = 'display: inline-block' )
      ),
    br(),
    HTML("<p> 2a. To check the selected genes for footprints of TFs or Pathways go to <a href='https://saezlab-funki-analysis-rnfy3j.streamlit.app/'> FUNKi </a> 
        from the <a href=' https://saezlab.org'> Saezlab </a and paste the genes into the gene submission field. 
        </p> "),
    HTML("<p> 2b. To enrich gene sets from various biological databases (GO, MSIG, KEGG etc.) go to <a href='https://maayanlab.cloud/Enrichr/'> Enrichr </a> 
        or <a href='https://maayanlab.cloud/enrichr-kg'> Enrichr-KG </a> from the <a href=' https://labs.icahn.mssm.edu/maayanlab/'> Mayan lab</a and paste the genes into the gene submission field. 
        </p> "),
    uiOutput("drugst_one_link"),
    p("4. Explore possible functional processes that your selected contrasts have in common and generate a hypothesis! "),
    hr(),  
    br(),
    br(),
    h4("Explore the signature in table form"),
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
              p("Choose the contrast to display functional footprinting results"),
              pickerInput(inputId = "select_tf",
                      label = "Select transcription factor(s)",
                      choices = sort(TFs),
                      multiple = T,
                      options = list(`live-search` = TRUE,
                                     size=6, `max-options` = 6),
                      selected = c("GATA4", "GATA3", "HIF1A", "HIF1B")
                      ),
              actionButton("reset_input_TF", "Reset TFs")
              ),
            mainPanel(
              tabsetPanel(
                type = "tabs",
                tabPanel("A. Animal models",
                         h3("Transcription factor activities"),
                         plotOutput("funcA_tf", width = "100%", height= "800px"),#, width = "100%", height = "600px"),
                         p("ribo, Ribo-seq (translational regulation); rna, RNA-seq (transcriptional regulation); 2d, two days; 2wk, two weeks;
                           swim, swimming (physiologic hypertrophy); tac, transverse-aortic-constriction (pahtologic hypertrophy)"),
                         br(),
                         br(),
                         hr()
                         ),
                tabPanel("B. Human Hypertrophic Cardiomyopathy",
                         ##Magnet
                         h3("Transcription factor activities - bulk HCM"),
                         plotOutput("funcB_tf_bulk", width = "100%", height = "500px"),
                         br(),
                         br(),
                         hr(),
                         h3("Transcription factor activities - single cell HCM"),
                         plotOutput("funcB_tf_sc", width = "100%", height = "500px"),
                         br(),
                         br(),
                         hr()
                         ),
                tabPanel("C+ D. Human Heart Failure and Fetal Program",
                         h3("Transcription factor activities"),
                        plotOutput("funcC_tf", width = "100%", height = "500px") %>%
                         withSpinner(),
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
                     includeMarkdown("inst/query_contrasts_sidebar.md"),
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
                                 selected = c("mm_TAC_RNA_2w") 
                     ),
                     
                     pickerInput(inputId = "select_contrast_hs_tf", 
                                 label = "B) Human HCM studies:",
                                 choices = sort(df_func%>% filter(cc =="B")%>% pull(condition)%>% unique()), 
                                 multiple = T,
                                 options = list(`live-search` = TRUE,
                                                size=10, `max-options` = 10),
                                 selected = c("hs_HCMvsNF_snRNA_CM", "hs_HCMvsNF_RNA") 
                     ),
                     
                     pickerInput(inputId = "select_contrast_hs2_tf", 
                                 label = "C+D) Human HF or fetal gene program",
                                 choices = df_func%>% filter(cc %in% c("C", "D"))%>% pull(condition)%>% unique(), 
                                 multiple = T,
                                 options = list(`live-search` = TRUE,
                                                size=10, `max-options` = 10),
                                 selected = c("hs_fetal_RNA") 
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
                     sliderInput("missing_prop_tf", "4. Select minimum number of contrasts to report a TF",
                                 min = 1, max = 10,
                                 value = 1, step= 1),
                     br(), 
                     strong("5. Compare contrasts!"), 
                     br(), 
                     actionButton(inputId = "submit_contrast_tf", label="Submit",
                                  icon=icon("paper-plane")) 
                   ),
                   mainPanel(
                     h3("Search for consistent genes"),
                     h4("Dysregulated genes (FDR-cutoff)"),
                     plotOutput("cq_hist_tf")%>%withSpinner(),
                     p("A. Barplot displays the overlap between the selected contrasts. Red indicates the number of genes for Signature generation."),
                     p("B. Barplot displays the number of genes of the intersect left that and their directionality of the intersection of all selected contrasts, i.e. how many genes are commonly up- or down- or inconsistently (up and down) regulated."),
                     p(""),
                     hr(),
                     br(),
                     br(),
                     h4("Full signature of dysregulated gene(s)"),
                     
                     plotOutput("hmap_top_tf")%>%withSpinner(),
                     hr(),
                     br(),
                     br(),
                     
                     h4("Top consistently dysregulated gene(s)"),
                     plotOutput("cq_top_tf")%>%withSpinner(),
                     hr(),
                     br(),
                     br()
                     )#main panel conserved tfs
                   ),# conserved TF panel
          tabPanel("Pathway activities", 
                   br(),
                   sidebarPanel(),
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
        title = "Enrichment analysis",
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
        title = "Acess Data",
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
          h4("How to acess full data"),
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
