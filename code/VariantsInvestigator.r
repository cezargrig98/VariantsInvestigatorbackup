library(shiny)
library(bslib)
library(shinyWidgets)
library(gt)
library(tibble)
library(DT)
library(DBI)
library(tidyverse)
# library(ggplot2)
library(ggbio)
library(plyranges)
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)
library(RSQLite)
library(Rsamtools)
library(trackViewer)
library(shinyFiles)



ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
 collapsible = TRUE,
      titlePanel("VariantsInvestigator"),
      div(img(
        src = "VariantInvestigator_logo_blue.png",
        height = 130,
        width = 200,
        style = "margin:1px 1px"
      ), ""),
      shinyFilesButton('files', label='File select', title='Please select an sql file', multiple=FALSE),
      
      sliderInput(
        inputId = "gnomad",
        "Maximum gnomAD allele frequency",
        min = 0,
        max = 1,
        step = 0.001,
        value = 1
      ),
      
      sliderInput(
        inputId = "SIFT",
        "Maximum SIFT value",
        min = 0,
        max = 1,
        step = 0.01,
        value = 1
      ),
      
      sliderInput(
        inputId = "PolyPhen",
        "Minimum PolyPhen value",
        min = 0,
        max = 1,
        step = 0.01,
        value = 0
      ),
      
      
      pickerInput(
        inputId = "impact",
        label = "Variant impact",
        choices = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
        selected = c("HIGH"),
        options = list(`actions-box` = TRUE),
        multiple = TRUE
      ),
      
      pickerInput(
        "HC",
        "High confidence LOF variants",
        choices = c("HC"),
        selected = c("HC"),
        multiple = TRUE
      ),
      
      pickerInput(
        inputId = "chr",
        label = "Chromosomes in analysis",
        choices = c(paste0("chr", as.character(c(1:22))),
                    "chrX",
                    "chrY",
                    "chrM"),
        selected = c(
          paste0(
            "chr", as.character(c(1:22))),
          "chrX",
          "chrY",
          "chrM"
        ),
        options = list(`actions-box` = TRUE),
        multiple = TRUE
      ),
      
      textInput("searchgene",
                "Filter by gene name",
                value = ""),
      
      pickerInput(
        inputId = "conseq",
        label = "Specific consequences (All consequences preselected)",
        choices = c(
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
          "sequence_variant"
        ),
        selected = c(
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
          "sequence_variant"
        ),
        options = c(`actions-box` = TRUE),
        multiple = TRUE
      ),
      
      # checkboxInput(inputId ="load", label = "Apply filters", value = FALSE)
      
      
    ),
    
    
    
    mainPanel(
      tabsetPanel(
        # card(
        tabPanel(
          #   card_header("Variants"),
          "Variants",
          
          
          
          
          dataTableOutput(outputId = "table")
          # full_screen = TRUE,
        ),
        
        
        tabPanel(
          # card_header(
          "Affected transcripts",
          # ),
          uiOutput("dropdownbutt"),
          
          
          plotOutput("genestbl")
        )
      )
    )
  )
)
server <- function(input, output, session) {
  output$myImage <- renderImage({
    list(
      src = "/www/VariantInvestigator_logo_blue.png",
      contentType = "image/png",
      width = 75,
      height = 50
    )
  }, deleteFile = FALSE)
  
  
  
  options(shiny.maxRequestSize=3000*1024^2)
  
  
  
  shinyFileChoose(input, 'files', root=c(root='.'), filetypes=c('sqlite', 'txt'))
  
  
  sqlfunc <- reactive({
    
    
    s_file <- parseFilePaths(root=c(root='.'), input$files)
    
    con <- dbConnect(RSQLite::SQLite(),
                     s_file$datapath
    )
    
    
    
    
    sqloutput <- dbSendQuery(con, paste("SELECT *
                          FROM variants
                          WHERE impact LIKE ?"),
                             params = list(c(paste(input$impact)))
    ) %>%
      dbFetch()  %>% separate_wider_delim("sift",
                                          delim = "(",
                                          names = c("sift_prediction","sift_values"),
                                          too_few = "align_start") %>%
      separate_wider_delim("polyphen",
                           delim = "(",
                           names = c("polyphen_prediction","polyphen_values"),
                           too_few = "align_start") %>% 
      dplyr::mutate(gnomad_af = as.numeric(gnomad_af))
    
    #risistemo nome della colonna SIFT per l'app finale
    sqloutput$sift_values <- gsub(")", "", sqloutput$sift_values) %>%
      as.numeric(sqloutput$sift_values)
    
    #stesso lavoro per polyphen
    sqloutput$polyphen_values <- gsub(")", "", sqloutput$polyphen_values) %>%
      as.numeric(sqloutput$polyphen_values)
    
    sqloutput  <- sqloutput %>%
      dplyr::filter(sift_values <= input$SIFT |
                      is.na(sift_values) == TRUE) %>%
      dplyr::filter(polyphen_values >= input$PolyPhen |
                      is.na(polyphen_values) == TRUE)
    
    
    sqloutput
    
    #   sqloutput <- sqloutput %>% dplyr::filter(grepl(paste(paste0("^", input$chr, "$"), collapse = "|"), chromosome)) %>%
    #     dplyr::filter(gnomad_af <= input$gnomad |
    #                     is.na(gnomad_af) == TRUE) %>%
    #     dplyr::filter(grepl(paste0("^", input$searchgene),
    #                         symbol,
    #                         ignore.case = TRUE)) %>%
    #     dplyr::filter(grepl(paste(input$conseq, collapse = "|"),
    #                         consequence)) %>%
    #     dplyr::filter(grepl(paste(input$HC, collapse = "|"), lof))
    #     
    # 
    
    
    
  })
  output$table <- renderDataTable(options = list(scrollX = TRUE,
                                                 width = "100%", height = "100%", responsive = TRUE),{
                                                   
                                                   #       if(input$load == TRUE) {
                                                   data_vars <- sqlfunc() %>% 
                                                     dplyr::filter(grepl(paste(paste0("^", input$chr, "$"), collapse = "|"), chromosome)) %>%
                                                     dplyr::filter(gnomad_af <= input$gnomad |
                                                                     is.na(gnomad_af) == TRUE) %>%
                                                     dplyr::filter(grepl(paste0("^", input$searchgene),
                                                                         symbol,
                                                                         ignore.case = TRUE)) %>%
                                                     dplyr::filter(grepl(paste(input$conseq, collapse = "|"),
                                                                         consequence)) %>%
                                                     dplyr::filter(grepl(paste(input$HC, collapse = "|"), lof)) %>%
                                                     dplyr::filter(sift_values <= input$SIFT |
                                                                     is.na(sift_values) == TRUE) %>%
                                                     dplyr::filter(polyphen_values >= input$PolyPhen |
                                                                     is.na(polyphen_values) == TRUE)
                                                   # }
                                                   
                                                   data_vars %>%
                                                     dplyr::mutate(symbol = paste0("<a target = '_blank' href = 'https://varsome.com/gene/hg38/", symbol, "'>",symbol,"</a>")) %>%
                                                     dplyr::mutate(symbol = map(symbol, gt::html))
                                                   
                                                   
                                                   
                                                   # gt() %>%
                                                   # fmt_number(decimals = 10,
                                                   #            drop_trailing_zeros = TRUE) %>%
                                                   # opt_interactive(
                                                   #   use_search = TRUE,
                                                   #   use_filters = TRUE,
                                                   #   use_resizers = TRUE,
                                                   #   use_highlight = TRUE,
                                                   #   use_compact_mode = TRUE,
                                                   #   use_text_wrapping = FALSE,
                                                   #   use_page_size_select = TRUE,
                                                   #   page_size_default = 10,
                                                   #   page_size_values = c(10, 25, 50, 100, 200)
                                                   # )
                                                   
                                                 })
  
  source("shinyfunctions.R")
  
  output$dropdownbutt <- renderUI({
    
    # if(input$load == TRUE) {
    data_vars <- sqlfunc() %>% 
      dplyr::filter(grepl(paste(paste0("^", input$chr, "$"), collapse = "|"), chromosome)) %>%
      dplyr::filter(gnomad_af <= input$gnomad |
                      is.na(gnomad_af) == TRUE) %>%
      dplyr::filter(grepl(paste0("^", input$searchgene),
                          symbol,
                          ignore.case = TRUE)) %>%
      dplyr::filter(grepl(paste(input$conseq, collapse = "|"),
                          consequence)) %>%
      dplyr::filter(grepl(paste(input$HC, collapse = "|"), lof)) %>% 
      dplyr::filter(sift_values <= input$SIFT |
                      is.na(sift_values) == TRUE) %>%
      dplyr::filter(polyphen_values >= input$PolyPhen |
                      is.na(polyphen_values) == TRUE)
    # }
    
    
    
    
    
    
    pickerInput(inputId = 'gene',
                label = 'Genes chosen',
                choices = c("",
                            unique(data_vars$symbol)
                )
    )
  })
  
  output$genestbl <- renderPlot({
    req(input$gene)
    
    sqlfunc() %>%
      dplyr::filter(symbol == input$gene) %>%
      gene_information_integration()
  })
  on.exit(dbDisconnect(con), add = TRUE)
}

shinyApp(ui, server)
