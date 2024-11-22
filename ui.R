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
                                     size=10, `max-options` = 6)
                      #selected = toupper(c("Nppb", "Nppa", "Mybpc3", "Col1a1", "Myh7", "Myh6" )) 
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
                     plotOutput("gene_expression_plots")%>%
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
                     h4("Explore posible associations of genetic variants with mouse phenotypes"),
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
                choices = unique(joint_contrast_df$contrast_id)[grep(pattern = "Mm|Rn", unique(joint_contrast_df$contrast_id))], 
                multiple = T,
                options = list(`live-search` = TRUE,
                               size=10, `max-options` = 10),
                selected = c("Mm_tac_ribo_2wk") 
    ),
    
    pickerInput(inputId = "select_contrast_hs", 
                label = "B) Human HCM studies:",
                choices = unique(joint_contrast_df$contrast_id)[grep(pattern = "HCM", unique(joint_contrast_df$contrast_id))], 
                multiple = T,
                options = list(`live-search` = TRUE,
                               size=10, `max-options` = 10),
                selected = c("Hs_singlecell_HCMvsNF_Cardiomyocyte", "Hs_bulk_HCMvsNF") 
    ),
    
    pickerInput(inputId = "select_contrast_hs2", 
                label = "C+D) Human HF or fetal gene program",
                choices = unique(joint_contrast_df$contrast_id)[!grepl(pattern = "HCM|Mm|Rn|DCM", unique(joint_contrast_df$contrast_id))], 
                multiple = T,
                options = list(`live-search` = TRUE,
                               size=10, `max-options` = 10),
                selected = c("Hs_fetal_Akat14") 
    ),
   
    actionButton("reset_input_contrasts", "Reset contrasts"),
    br(),
    #br(),
    hr(),
    pickerInput(inputId = "select_alpha", 
                label = "2. Select alpha for FDR cut off:",
                choices = c(0.0001, 0.001, 0.01, 0.05, 0.1), 
                multiple = F,
                options = list(`live-search` = TRUE,
                               size=10, `max-options` = 10),
                selected = 0.05
                ),
    
    sliderInput("cut_off_genes", "3. Select number of genes to plot:",
                min = 1, max = 40,
                value = 10, step= 1),
    hr(),
    sliderInput("missing_prop", "4. If you want to allow for NAs, select number of contrasts with NAs allowed",
                min = 0, max = 10,
                value = 0, step= 1),
    br(), 
    strong("4. Compare contrasts!"), 
    br(), 
    actionButton("submit_contrast", label="Submit",
                 icon=icon("paper-plane")) 
  ),
  mainPanel(
    h3("Search for consistent genes"),
    h4("Dysregulated genes (FDR-cutoff)"),
    plotOutput("cq_hist")%>%withSpinner(),
    p("A. Venn Diagram displays intersections of gene sets from selected contrasts at chosen FDR cut off"),
    p("B. Barplot displays number of genes regarding their directionality of the intersection of all selected contrasts, i.e. how many genes are commonly up- or down- or inconsistently (up and down) regulated."),
    p(""),
    hr(),
    br(),
    br(),
    
    h4("Top consistently dysregulated gene(s)"),
    plotOutput("cq_top")%>%withSpinner(),
    hr(),
    br(),
    br(),
    
    h4("Top consistently dysregulated gene(s)"),
    plotOutput("hmap_top")%>%withSpinner(),
    hr(),
    br(),
    br(),
    
    h4("Table of consistently dysregulated gene(s)"),
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
    
    h3("Characterize conistent genes"),
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
    HTML("<p> 2b. To enrich gene sets from various biological databases (GO, MSIG, KEGG etc.) go to <a href='https://maayanlab.cloud/enrichr-kg'> Enrichr </a> 
        from the <a href=' https://labs.icahn.mssm.edu/maayanlab/'> Mayan lab</a and paste the genes into the gene submission field. 
        </p> "),
    HTML("<p> 2c. To explore possible gene network based drug interactions, tissue expression and related disorders go to <a href='https://drugst.one/standalone'> Drugst.one </a> f
      from th the <a href='https://www.cosy.bio/'> Baumbachlab </a> and paste the genes into the network input field.
        </p> "),
    p("3. Explore possible functional processes that your selected contrasts have in common and generate a hypothesis! "),
    br(),
    br(),
    hr()
     
  )
),

## old, used drugst.one intersect
# tabPanel(
#   title = "Test drugst.one",
#   icon = icon("chart-line"),
#   HTML( "<drugst-one id='drugstone-component-id'></drugst-one>")
# ),
      
      #### Functional analysis ####
tabPanel(
        title = "Functional Genomics",
        icon = icon("chart-line"),
        sidebarPanel(
          includeMarkdown("inst/functional_analysis_sidebar.md"),


          p("Choose the contrast to display functional footprinting results"),
        
          pickerInput(inputId = "select_tf",
                      label = "Select transcription factor(s)",
                      choices = sort(TFs),
                      multiple = T,
                      options = list(`live-search` = TRUE,
                                     size=6, `max-options` = 6)
                      #selected = c("GATA4", "GATA3", "HIF1A", "HIF1B")
          ),
          actionButton("reset_input_TF", "Reset TFs")

        ),
        mainPanel(
          tabsetPanel(
            type = "tabs",
            tabPanel("A. Animal models",
                     h3("Transcription factor activities"),
                     plotOutput("funcA_tf"),#, width = "100%", height = "600px"),
                     p("ribo, Ribo-seq (translational regulation); rna, RNA-seq (transcriptional regulation); 2d, two days; 2wk, two weeks;
               swim, swimming (physiologic hypertrophy); tac, transverse-aortic-constriction (pahtologic hypertrophy)"),
                     br(),
                     br(),
                     hr(),
                     h3("Pathway activities"),
                     plotOutput("funcA_pw"),#, width = "100%", height = "600px"),
                     br(),
                     br(),
                     hr(),
                      tabsetPanel(
                        type = "tabs",
                        tabPanel("TF-activities",
                                 DT::dataTableOutput("funcA_tb_tf")),
                         tabPanel("Pathway-activities",
                                DT::dataTableOutput("funcA_tb_pw"))
                      )

            ),
            tabPanel("B. Human Hypertrophic Cardiomyopathy",
                     ##Magnet
                     h3("Transcription factor activities - bulk HCM"),
                     plotOutput("funcB_tf_bulk", width = "100%", height = "500px"),
                     br(),
                     
                     h3("Pathway activities - bulk HCM"),
                     plotOutput("funcB_pw_bulk", width = "100%", height = "500px"),
                     
                     br(),
                     br(),
                     hr(),
                     h3("Transcription factor activities - single cell HCM"),
                     plotOutput("funcB_tf_sc", width = "100%", height = "500px"),
                     br(),
                     br(),
                     hr(),
                     h3("Pathway activities - single cell HCM"),
                     plotOutput("funcB_pw_sc", width = "100%", height = "500px"),
                     p("HCM, hypertrophic cardiomyopathy; DCM, dilated cardiomyopathy; NF, non-failing"),
                     br(),
                     br(),
                     hr(),
                     tabsetPanel(
                       type = "tabs",
                       tabPanel("TF-activities bulk",
                                DT::dataTableOutput("funcB_tf_tb_bulk")),
                       tabPanel("TF-activities single cell ",
                                DT::dataTableOutput("funcB_tf_tb_sc")),
                       tabPanel("Pathway-activities bulk",
                                DT::dataTableOutput("funcB_pw_bulk_tb")),
                       tabPanel("Pathway-activities single cell",
                                DT::dataTableOutput("funcB_pw_sc_tb"))
                     )
            ),
            tabPanel("C. Human Heart Failure",
                     h3("Transcription factor activities"),
                    plotOutput("funcC_tf", width = "100%", height = "500px") %>%
                     withSpinner(),
                    br(),
                    br(),
                    hr(),
                    h3("Pathway activities"),
                     plotOutput("funcC_pw", height= "250px"),
                    tabsetPanel(
                      type = "tabs",
                      tabPanel("TF-activities",
                               DT::dataTableOutput("funcC_tf_tb")),
                      tabPanel("Pathway-activities",
                               DT::dataTableOutput("funcC_pw_tb"))
                    )
                     ),
            tabPanel("D. Fetal gene program",
                     h3("Transcriptionfactor activities"),
                     plotOutput("funcD_tf", height= "250px"),
                     br(),
                     br(),
                     hr(),
                     h3("Pathway activities"),
                     plotOutput("funcD_pw", height= "250px"),
                     p("rna_fetal1, fetal study 1 (Spurell et al, 2022); rna_fetal2, fetal study 2 (Akat et al, 2014)"),
                     tabsetPanel(
                       type = "tabs",
                       tabPanel("TF-activities",
                                DT::dataTableOutput("funcD_tf_tb")),
                       tabPanel("Pathway-activities",
                                DT::dataTableOutput("funcD_pw_tb"))
                     )
                     )
            )
          )
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
        title = "Download Center",
        icon = icon("table"),
        sidebarPanel(
          includeMarkdown("inst/meta_analysis_results_sidebar.md")
        ),
        mainPanel(
          h5("Download gene expression contrast data"),
          h6("Select your desired contrast, explore and sort data in table format or download"),
          tabsetPanel(
            type = "tabs",
            tabPanel("A. Animal Models", DT::dataTableOutput("mouse_hypertrophyDT"),
                     br(), p("Download full data from zenodo, LINK")),
            tabPanel("B. Human HCM", DT::dataTableOutput("human_HCMDT"),
                     br(), p("Download full data from zenodo, LINK")),
           tabPanel("C. Human HF",  DT::dataTableOutput("human_HF_bulk_summDT"),
                     br(), p("Download full data from zenodo, LINK")),
            tabPanel("D. Fetal gene program", DT::dataTableOutput("human_fetalDT"),
                     br(), p("Download full data from zenodo, LINK"))
          ),
          br(), 
          hr(),
          
        )
      )
      
      #### Footer ####
      #footer = column(12, align="center", "CHEERIO-App")
    ) # close navbarPage
  ) # close fluidPage
}
