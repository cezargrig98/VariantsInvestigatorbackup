library(shiny)
library(bslib)
library(shinyWidgets)
library(gt)
library(tibble)
library(DT)
library(DBI)
library(tidyverse)
library(plyranges)
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(RSQLite)
library(Rsamtools)
library(trackViewer)
library(shinyFiles)

#theme information
light <- bs_theme(version = 5, bootswatch = "zephyr", primary = "rgb(8, 14, 116, 255)", 
                  info = "rgb(8, 14, 116, 255)", warning = "rgb(165, 166, 246, 255)")
dark <- bs_theme(version = 5, preset = "zephyr", bg = "rgb(0, 0, 0)", fg = "rgb(255, 255, 255)", 
                 primary = "rgb(165, 166, 246, 255)", secondary = "#287BB5", 
                 success = "rgb(165, 166, 246, 255)", info = "rgb(165, 166, 246, 255)", 
                 warning = "rgb(165, 166, 246, 255)", danger = "rgb(165, 166, 246, 255)")
#ui part
ui <- page_sidebar(
  sidebar = sidebar(
    div(img(
      src = "VariantInvestigator_logo_blue.png",
      height = 130,
      width = 200,
      style = "margin:1px 1px"
    ), ""),
    
    #standard theme
    theme = light,
    
    #night mode switch
    materialSwitch
    (
      "dark_mode",
      "Dark mode",
      status = "primary"
    ),
    
    #input file
    shinyFilesButton
    (
      'files', 
      label='Select file', 
      title='Please select an SQL file', 
      multiple=FALSE, 
      buttonType = "primary"
      ),
    
    #lof filter
    div
    (
      pickerInput
      (
        "HC",
        "LoF variants",
        choices = c("HC"),
        selected = c(""),
        multiple = TRUE
      ),
      shiny::tags$small(style = "color: gray; font-size: 12px;", "High-confidence loss-of-function variants")
    ),
    
    #impact filter
    pickerInput
    (
      inputId = "impact",
      label = "Variant impact",
      choices = c("HIGH", "MODERATE", "LOW", "MODIFIER", "NA"),
      selected = c("HIGH", "MODERATE", "NA"),
      options = list(`action-box` = TRUE),
      multiple = TRUE
    ), 
    
    #clinical significance filter
    pickerInput
    (
      inputId = "clinsign",
                label = "ClinVar classifications",
                choices = c("benign",
                            "uncertain_significance",
                            "pathogenic",
                            "drug_response",
                            "association",
                            "affects",
                            "other",
                            "conflicting",
                            "NA"),
                selected = c("pathogenic", "conflicting", "NA"),
                options = c(`actions-box` = TRUE),
                multiple = TRUE
    ),
    
    #consequence filter
    pickerInput(
      inputId = "conseq",
      label = "Variant consequences",
      choices = 
        c(
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_amplification",
        "feature_elongation",
        "feature_truncation",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        "splice_donor_5th_base_variant",
        "splice_region_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
        "incomplete_terminal_codon_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "regulatory_region_variant",
        "intergenic_variant",
        "sequence_variant",
        "NA"
      ),
      selected = 
        c(
        "transcript_ablation",
        "splice_acceptor_variant",
        "splice_donor_variant",
        "stop_gained",
        "frameshift_variant",
        "stop_lost",
        "start_lost",
        "transcript_amplification",
        "feature_elongation",
        "feature_truncation",
        "inframe_insertion",
        "inframe_deletion",
        "missense_variant",
        "protein_altering_variant",
        "splice_donor_5th_base_variant",
        "splice_region_variant",
        "splice_donor_region_variant",
        "splice_polypyrimidine_tract_variant",
        "incomplete_terminal_codon_variant",
        "start_retained_variant",
        "stop_retained_variant",
        "synonymous_variant",
        "coding_sequence_variant",
        "mature_miRNA_variant",
        "5_prime_UTR_variant",
        "3_prime_UTR_variant",
        "non_coding_transcript_exon_variant",
        "intron_variant",
        "NMD_transcript_variant",
        "non_coding_transcript_variant",
        "coding_transcript_variant",
        "upstream_gene_variant",
        "downstream_gene_variant",
        "TFBS_ablation",
        "TFBS_amplification",
        "TF_binding_site_variant",
        "regulatory_region_ablation",
        "regulatory_region_amplification",
        "regulatory_region_variant",
        "intergenic_variant",
        "sequence_variant",
        "NA"
      ),
      options = c(`actions-box` = TRUE),
      multiple = TRUE
    ),
    
    #chromosome choice filter
    pickerInput
    (
      inputId = "chr",
      label = "Chromosomes in analysis",
      choices = c(paste0("chr", as.character(c(1:22))),
                  "chrX",
                  "chrY",
                  "chrM"),
      selected = 
        c(
          paste0(
            "chr", as.character(c(1:22))),
          "chrX",
          "chrY",
          "chrM"
        ),
      options = list(`actions-box` = TRUE),
      multiple = TRUE
    ),
    
    #gnomad filter
    sliderInput
    (
      inputId = "gnomad",
      "GnomAD allele frequency",
      min = 0,
      max = 1,
      step = 0.05,
      value = 0.5
    ),
    
    #allele depth filter
    sliderInput
    (
      inputId = "dp",
      "Read depth", #Minimum and maximum depth read
      min = 0, max = 200,
      value = c(0, 200)
    ),
    
    #sift filter
    sliderInput
    (
      inputId = "SIFT",
      "SIFT value", #Maximum SIFT value
      min = 0,
      max = 1,
      step = 0.01,
      value = 1
    ),
    
    #polyphen filter
    sliderInput
    (
      inputId = "PolyPhen",
      "PolyPhen value", #Minimum 
      min = 0,
      max = 1,
      step = 0.01,
      value = 0
    ),
     
   #gene search filter
    div(
      textInput
      (
        "searchgene", 
        "Gene name(s)", 
        value = ""
        ),
      shiny::tags$small(style = "color: gray; font-size: 12px;", "Use `|` separator for multiple genes")
    )
   
  ),
  
  #tabs
  navset_tab(
    
    #variant table tab
    nav_panel(
      "Variants",
      DTOutput(outputId = "table"),
      fillable = FALSE
    ),
    
    #genome view tab
    nav_panel(
      "Affected transcripts",
      uiOutput("dropdownbutt"),


      plotOutput("genestbl")
    ),
    
    #number of mutations per gene tab
    nav_panel(
      "Mutations on each gene",
      DTOutput(outputId = "symbol_num")
      
    )
  )
)

#server processing part
server <- function(input, output, session) {
  

  
  observe(session$setCurrentTheme(
    if (isTRUE(input$dark_mode)) dark else light
  ))

  #app logo
  output$myImage <- renderImage({
    list(
      src = "/www/VariantInvestigator_logo_blue.png",
      contentType = "image/png",
      width = 75,
      height = 50
    )
  }, deleteFile = FALSE)
  
  #allowing big files in shiny
  options(shiny.maxRequestSize=3000*1024^2)
  
  #file input server side processing
  shinyFileChoose(input, 'files', root=c(root='.'), filetypes=c('sqlite', 'txt'))
  
  #reactive function for database creation
  sqlfunc <- reactive({
    
  #error message for no file input
    s_file <- parseFilePaths(
      root=c(root='.'), input$files)
     validate(
       need(s_file$datapath != "", message = "Please select a data set"
       )
     ) 
     
    #sql database connection
    con <- dbConnect(RSQLite::SQLite(),
                     s_file$datapath
    )
    
    #first filter and data loading in memory
    sqloutput <-  dbSendQuery(con, paste("SELECT *
                          FROM variants
                          WHERE impact LIKE ?"),
                              params = list(c(paste(input$impact)))
    ) %>%
      dbFetch()
    
    #genotype table fetching
    Genotype <- dbGetQuery(con, "SELECT AD, DP, GT, start
                     FROM genotypes")
    
    #genotype and variants merging in single table
    sqloutput <- merge(sqloutput, Genotype, by.x="start", by.y="start", all.x = TRUE) %>% 
      separate_wider_delim("sift",
                           delim = "(",
                           names = c("sift_prediction","sift_values"),
                           too_few = "align_start") %>%
      separate_wider_delim("polyphen",
                           delim = "(",
                           names = c("polyphen_prediction","polyphen_values"),
                           too_few = "align_start") %>%
      dplyr::mutate(gnomad_af = as.numeric(gnomad_af)) %>%
      relocate("ad", .before = "allele") %>%
      relocate("dp", .before = "ad") %>%
      relocate("gt", .before = "allele") %>% 
      relocate("start", .before = "end") %>% 
      dplyr::select(-"ac")
    
    #conversion of certain columns to facilitate browsing
    sqloutput$dp <- as.numeric(sqloutput$dp)
    sqloutput$ad <- gsub(",", ", ", sqloutput$ad)
    
    
    #conversion of sift column to numeric
    sqloutput$sift_values <- gsub(")", "", sqloutput$sift_values) %>%
      as.numeric(sqloutput$sift_values)
    
    #same for polyphen
    sqloutput$polyphen_values <- gsub(")", "", sqloutput$polyphen_values) %>%
      as.numeric(sqloutput$polyphen_values)
    
    #converting '.' values to "NA" to avoid annoyance in different representation of missing values
    sqloutput$clin_sig <- gsub("\\.", "NA", sqloutput$clin_sig)
    
    #applying filtersd
      sqloutput <- sqloutput %>%
        dplyr::filter(grepl(paste(input$HC, collapse = "|"),
                            lof)) %>% 
        dplyr::filter(sift_values <= input$SIFT |
                        is.na(sift_values) == TRUE) %>%
        dplyr::filter(polyphen_values >= input$PolyPhen |
                        is.na(polyphen_values) == TRUE) %>%
        dplyr::filter(dp >= input$dp[1] & dp <= input$dp[2] | is.na(dp) == TRUE) %>% 
        dplyr::filter(grepl(paste(paste0("^", input$chr, "$"), collapse = "|"), chromosome)) %>%
        dplyr::filter(gnomad_af <= input$gnomad |
                        is.na(gnomad_af) == TRUE) %>%
        dplyr::filter(grepl(paste0("^", input$searchgene),
                            symbol,
                            ignore.case = TRUE)) %>%
        dplyr::filter(grepl(paste(input$conseq, collapse = "|"),
                            consequence)) %>%
        dplyr::filter(grepl(paste(c(input$clinsign), collapse = "|"),
                            clin_sig))
      
      #result of function
      sqloutput

    
  })
  
  #variants table output
  output$table <- renderDT({
    
    #correct display of variants table
    data_vars <- sqlfunc() %>% 
      dplyr::select(-"gene") %>% 
      dplyr::mutate(symbol = paste0("<a target = '_blank' href = 'https://varsome.com/gene/hg38/", symbol, "'>",symbol,"</a>")) %>%
      dplyr::mutate(symbol = map(symbol, gt::html)) %>% 
      relocate(symbol) %>% 
      relocate(consequence, .before = id) %>% 
      relocate(clin_sig, .before = chromosome) %>%
      rowwise() %>%
      mutate(ref = case_when(nchar(ref) > 10 ~
                               paste(str_sub(ref, 1, 10), "..."),
                             nchar(ref) <= 10 ~ ref)) %>%
      mutate(alt = case_when(nchar(alt) > 10 ~
                               paste(str_sub(alt, 1, 10), "..."),
                             nchar(alt) <= 10 ~ alt)) %>% 
      mutate(clin_sig = case_when(nchar(clin_sig) > 20 ~
                               paste(str_sub(clin_sig, 1, 20), "..."),
                             nchar(clin_sig) <= 20 ~ clin_sig))
    
    DT::datatable(
      data = data_vars,
      class = 'cell-border stripe', 
      plugin = "ellipsis",
      filter = 'top', selection = 'none',
      height = "100%", width = "100px",
      
    )
  }, server = TRUE)
  
  #graphical tab function
  source("shinyfunctions.R")
  
  #gene selector
  output$dropdownbutt <- renderUI({
    
   #ordering genes in alphabetical order for convenience
    data_vars <- sqlfunc() %>%
      dplyr::select(symbol) %>% 
      arrange(symbol)
    
    #selector of gene for graphical representation
    pickerInput(inputId = 'gene',
                label = 'Genes chosen',
                choices = c("",
                            unique(data_vars$symbol)
                )
    )
  })
  
  output$symbol_num <- renderDT({
    
    #number of variants per gene table
    data_vars <- sqlfunc() %>%
      dplyr::select(symbol) %>% 
      dplyr::group_by(symbol) %>% 
      dplyr::summarise(total_vars = dplyr::n(),
                       .groups = 'drop') %>% 
      dplyr::arrange(desc(total_vars)) 

    
    DT::datatable(
      data = data_vars,
      class = 'cell-border stripe', 
      plugin = "ellipsis",
      filter = 'top', selection = 'none',
      height = "100%", width = "100px"
    )
    
  }, server=TRUE)

  #gene rendering functions
  output$genestbl <- renderPlot({
    req(input$gene)

    sqlfunc() %>%
      dplyr::filter(symbol == input$gene) %>%
      gene_information_integration()
  })
  
}

shinyApp(ui, server)
