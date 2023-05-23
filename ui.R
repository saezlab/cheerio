# UI
source("sub/global.R")
source("sub/helper.R")
ui = function(request) {
  fluidPage(
    useShinyjs(),
    #tags$head(includeScript("google-analytics.js")),
    navbarPage(
      id = "menu", 
      title = div(img(src="logo_saezlab.png", width="25", height="25"),
                  "CHROMA"),
      windowTitle = "Murine Cardiac Hypertrophy",
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
                      choices = sort(unique(contrasts$MgiSymbol)), 
                      multiple = T,
                      options = list(`live-search` = TRUE,
                                     size=10, `max-options` = 10),
                      selected = c("Nppb", "Nppa", "Postn", "Col1a1", "Myh7", "Myh6" ) 
                      ),
          
          actionButton("reset_input", "Reset genes"),
          
          checkboxGroupInput(inputId= "contrasts", label= "Select contrasts", 
                             selected = c("Murine_Hypertrophy"),
                             inline = FALSE, width = NULL, 
                             choiceNames = c("Murine Hypertrophy" , "Human heart failure", "Human fetal" ),
                             choiceValues = c("Murine_Hypertrophy" , "Human_HF", "Fetal" ))
        ),
        mainPanel(
          h3("Regulation in Murine Cardiac Hypertrophy models"),
          h5("Gene expression regulation"),
          plotOutput("gene_expression_plots"),#, width = "100%", height = "600px"),
          h5("Ribo-seq and RNA-seq correlation"),
          plotOutput("cardiac_hyper_corr"),
          h3("Regulation in Human Heart Failure studies"),
          h5("Regulation in single bulk studies"),
          h6("HF bulk transcripomic studies"),
          plotOutput("HFgene_regulation_boxplot", width = "100%", height = "500px") %>%
            withSpinner(),
          h5("Distribution of mean t-values"),
          h6("Mean t-values distribution of all bulk studies"),
          plotlyOutput("mean_t_dist", width = "100%", height = "250px") %>%
            withSpinner(),
          h5("Ranking of queried genes"),
          h6("Consensus ranking, the higher the rank the more consistently is the gene regulated in human HF"),
          plotlyOutput("rank_position", width = "100%", height = "125px") %>%
            withSpinner(),
        
          h5("Single cell gene expression"),
          h6("Celltype regulation compared between non failing, DCM and HCM"),
          plotOutput("HF_single", width = "100%") %>%
            withSpinner(),
          
          h3("Regulation in Fetal Hearts"),
          plotOutput("fetal_gene_expression_plots"),
        
          hr()
         
        )
      ),
      
      

# Query contrasts -----------------------------------------------------------------------------
tabPanel(
  title = "Query contrasts",
  icon = icon("search"),
  sidebarPanel(
    includeMarkdown("inst/query_contrasts_sidebar.md"),
    pickerInput(inputId = "select_contrast", 
                label = "Select contrasts",
                choices = sort(unique(joint_contrast_df$contrast_id)), 
                multiple = T,
                options = list(`live-search` = TRUE,
                               size=10, `max-options` = 10),
                selected = c("tac_ribo_2wk", "fetal_rna_fetal1", "HCMvsNF_Fibroblast") 
    ),
    
    actionButton("reset_input_contrasts", "Reset contrasts"),
    
    sliderInput("cut_off_genes", "Cut off top genes:",
                min = 1, max = 30,
                value = 10, step= 1),
    
    pickerInput(inputId = "select_alpha", 
                label = "Select alpha for FDR cut off",
                choices = c(0.0001, 0.001, 0.01, 0.05, 0.1), 
                multiple = F,
                options = list(`live-search` = TRUE,
                               size=10, `max-options` = 10),
                selected = 0.05
    ),
    
    p("Contrast comparison performed clicking the submit button."),
    actionButton("submit_contrast", label="Submit",
                 icon=icon("send")) 
  ),
  mainPanel(
    h3("Search for consistent genes"),
    h5("Hist"),
    plotOutput("cq_hist"),
    
    h5("Top hits"),
    plotOutput("cq_top"),
    
    h5("table"),
    DT::dataTableOutput("cq_table")
    )
    
  
),

      
      #### Functional analysis ####
      tabPanel(
        title = "Functional analysis",
        icon = icon("chart-line"),
        sidebarPanel(
          includeMarkdown("inst/functional_analysis_sidebar.md"),
          
          h3(),
          p("Choose the contrast to display functional footprinting results"),
          # radioButtons("select_contrast_func2", "Select contrast", 
          #              choices = c("murine_hypertrophy", 
          #                          "human_HF", 
          #                          "human_fetal")),
          selectizeInput("select_contrast_func", "Select contrast", 
                         multiple = F, 
                         choices = c("murine_hypertrophy", 
                                     "human_HF",
                                     "human_HF_sc",
                                     "human_fetal")),
          h2(),
          pickerInput(inputId = "select_tf", 
                      label = "Select transcription factor(s)",
                      choices = sort(unique(tf_hypertrophy$source)), 
                      multiple = T,
                      options = list(`live-search` = TRUE,
                                     size=6, `max-options` = 6),
                      selected = c("Gata4", "Gata3", "Hif1a")
          ),
          actionButton("reset_input_TF", "Reset TFs")
        
        ),
        mainPanel(
          h3("Estimated TF activity"),
          h5(""),
          plotOutput("tf_hypertrophy_plot",width = "100%", height = "900px"),
          
          h3("Estimated Pathway activity"),
          plotOutput("progeny_hypertropy_plot"),

          h3("Access, query and download raw data"),
          h5(""),
          tabsetPanel(
            type = "tabs",
            tabPanel("DoRothEA",
                     DT::dataTableOutput("dorothea_table_hypertrophy")),
            tabPanel("PROGENy", DT::dataTableOutput("progeny_table"))
            
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
          checkboxGroupInput(inputId= "contrasts_gsea", label= "Select contrasts to perform enrichment on", 
                             selected = c("Murine_Hypertrophy"),
                             inline = FALSE, width = NULL, 
                             choiceNames = c("Murine Hypertrophy" , "Human heart failure bulk", 
                                             "Human heart failure sc", "Human fetal" ),
                             choiceValues = c("murine_hypertrophy" , "human_HF","human_HF_sc",  "fetal" )),
          
          p("Choose whether you would like to test your gene set agains a 
             directed or undirected signature."),
          radioButtons("signature_source", "Signature", 
                       choices = c("directed", "undirected")),
          p("GSEA will be performed upon clicking the submit button."),
          actionButton("submit", label="Submit",
                       icon=icon("send")) 
          ),
        mainPanel(
          h4("GSEA result"),
          DT::dataTableOutput("gsea_res_table"),
          hr(),
          h4("GSEA plots"),
          plotOutput("gsea_res_plots", width = "100%", height = "600px") %>%
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
            tabPanel("Murine Hypertrophy", DT::dataTableOutput("mouse_hypertrophyDT")),
            tabPanel("Human HF Bulk individual",  DT::dataTableOutput("human_HF_bulk_indDT")),
            tabPanel("Human HF Bulk summary",  DT::dataTableOutput("human_HF_bulk_summDT")),
            tabPanel("Human HF single cell", DT::dataTableOutput("human_scDT")),
            tabPanel("Human fetal bulk", DT::dataTableOutput("human_fetalDT"))
          ),
          h6("*Human HF Bulk individual: Contains single study results from ReHeaT"),
          h6("*Human HF Bulk summary: Contains consensus ranking (summary of single studies) results from ReHeaT")
        )
      ),
      
      
      #### Study overview ####
      # tabPanel(
      #   title = "Study overview",
      #   icon = icon("database"),
      #   sidebarPanel(
      #     includeMarkdown("inst/overview_sidebar.md")
      #   ),
      #   mainPanel(
      #     DT::dataTableOutput("overview"),
      #     br(),
      #     includeMarkdown("inst/overview.md")
      #   )
      # ),
      
      #### Footer ####
      footer = column(12, align="center", "CH-App 2020")
    ) # close navbarPage
  ) # close fluidPage
}
