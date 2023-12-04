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



ui <- page_sidebar(
  title = "VariantsInvestigator",
  sidebar = sidebar(
    title = "Filters",
    sliderInput(
      inputId = "gnomad",
      "Maximum gnomAD allele frequency",
      min = 0,
      max = 1,
      step = 0.001,
      value = 1
    ),
    # sliderInput(
    #   inputId = "geno_qual",
    #   "minimum genotype quality",
    #   min = 0,
    #   max = 99,
    #   step = 1,
    #   value = 0
    # ),
    # sliderInput(
    #   inputId = "SIFT",
    #   "Maximum SIFT value",
    #   min = 0,
    #   max = 1,
    #   step = 0.01,
    #   value = 1
    # ),
    # sliderInput(
    #   inputId = "PolyPhen",
    #   "Minimum PolyPhen value",
    #   min = 0,
    #   max = 1,
    #   step = 0.01,
    #   value = 0
    # ),

    pickerInput(
      inputId = "impact",
      label = "Variant impact",
      choices = c("HIGH", "MODERATE", "LOW", "MODIFIER"),
      selected = c("HIGH"),
      options = list(`actions-box` = TRUE),
      multiple = TRUE
    ),
    
    pickerInput(
      inputId = "chr",
      label = "Chromosomes in analysis",
      choices = c(paste0("chr", as.character(c(1:22))),
                  "chrX",
                  "chrY",
                  "chrM"),
      selected = c(paste0("chr", as.character(c(1:22))),
                   "chrX",
                   "chrY",
                   "chrM"),
      options = list(`actions-box` = TRUE),
      multiple = TRUE
    ),
    # 
    textInput("searchgene",
              "Filter by gene name",
              value = ""),
    # 
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
    # 
    pickerInput(
      "HC",
      "High confidence LOF variants",
      choices = c("HC"),
      selected = c("HC"),
      multiple = TRUE
    )
    
  ),
  layout_columns(
    card(
      card_header("Variants"),
      gt_output(outputId = "table"),
      full_screen = TRUE
    )
    ),
  layout_columns(
    card(
      card_header("Affected genes"),
      full_screen = TRUE,
          
          dropdownButton(
            

            
            pickerInput(inputId = 'chr_plot',
                        label = 'Chromosomes displayed',
                        choices = c(paste0("chr", as.character(c(1:22))),
                                    "chrX",
                                    "chrY",
                                    "chrM"),
                        selected = c("chr1"),
                        options = c(`actions-box` = TRUE),
                        multiple = TRUE),
            
            circle = TRUE,
            icon = icon("gear"), width = "50px"
          ),
          
          plotOutput("chroms")
      ),
      width = "400px",
  card(
    card_header("Position of the variant on the gene"),
    full_screen = TRUE,
    plotOutput("gene")
      ),
  width = "400px"
  )
)

server <- function(input, output, session) {
  # output$myImage <- renderImage({
  #   list(
  #     src = "/home/shared_projects/TESI/tesi_cezar_grigorean/code/www/VariantInvestigator_logo_blue.png",
  #     contentType = "image/png",
  #     width = 150,
  #     height = 100
  #   )
  # }, deleteFile = FALSE)

 
  
    output$table <- render_gt({
    req(input$impact)
     con <- dbConnect(RSQLite::SQLite(),
              "/home/shared_projects/TESI/tesi_cezar_grigorean/data/test.sqlite"
              ) 
     on.exit(dbDisconnect(con), add = TRUE) 
    dbSendQuery(con, paste("SELECT *
                          FROM variants
                          WHERE impact LIKE ?"),
                params = list(c(paste(input$impact)))
                        ) %>%
        dbFetch() %>%
      dplyr::filter(grepl(paste(paste0("^", input$chr, "$"), collapse = "|"), chromosome)) %>%
      
      dplyr::filter(gnomad_af <= input$gnomad |
                      is.na(gnomad_af) == TRUE) %>%
      dplyr::filter(grepl(paste0("^", input$searchgene),
                          symbol,
                          ignore.case = TRUE)) %>%
      dplyr::filter(grepl(paste(input$conseq, collapse = "|"),
                          consequence)) %>%
      dplyr::filter(grepl(paste(input$HC, collapse = "|"), lof)) %>%
      gt() %>%
      fmt_number(decimals = 10,
                 drop_trailing_zeros = TRUE) %>%
      opt_interactive(
        use_search = TRUE,
        use_filters = TRUE,
        use_resizers = TRUE,
        use_highlight = TRUE,
        use_compact_mode = TRUE,
        use_text_wrapping = FALSE,
        use_page_size_select = TRUE,
        page_size_default = 10,
        page_size_values = c(10, 25, 50, 100, 200)
      )

  })
    
    chr_information_integration <- function(x){
      
      
      # grch37dict = read_tsv("/home/shared_projects/REFS/Homo_sapiens/GATK/GRCh37/Sequence/WholeGenomeFasta/human_g1k_v37_decoy.dict",
      #                       skip = 1, col_names = c("class", "sequence", "length", "m5", "as", "url", "species"))
      
      
      grch38dict = read_tsv("/home/shared_projects/REFS/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict",
                            skip = 1, col_names = c("class", "sequence", "length", "m5", "as", "url", "species"))
      # 
      # grch37dict$sequence <- gsub("SN:", "", grch37dict$sequence)
      # grch37dict$length <- gsub("LN:", "", grch37dict$length)
      
      grch38dict$sequence <- gsub("SN:", "", grch38dict$sequence)
      grch38dict$length <- gsub("LN:", "", grch38dict$length)
      
      # newSeqInfoData37 <- Seqinfo(
      #   seqnames = grch37dict$sequence,
      #   seqlengths = as.numeric(grch37dict$length),
      #   isCircular = NA,
      #   genome = "hg37"
      # )
      
      
      # txdb <- makeTxDbFromGFF("/home/shared_projects/REFS/Homo_sapiens/NCBI/GRCh38/Annotation/Genes/genes.gtf")
      
      
      # geneGranges <- transcripts(txdb)
      
      
      
      
      newSeqInfoData38 <- Seqinfo(
        seqnames = grch38dict$sequence,
        seqlengths = as.numeric(grch38dict$length),
        isCircular = NA,
        genome = "hg38"
      )
      # 
      # newSeqInfoData37_noextra <- newSeqInfoData37[
      #   c(
      #     as.character(c(1:22)),
      #     "X",
      #     "Y",
      #     "MT"
      #   )
      # ]
      
      newSeqInfoData38_noextra <- newSeqInfoData38[
        c(
          paste0("chr", as.character(c(1:22))),
          "chrX",
          "chrY",
          "chrM"
        )
      ]
      
      integration <- GRanges(
        x$chromosome,
        IRanges(
          start = x$start,
          end = x$end
        )
      )
      
      integration_noextra <- integration[
        seqnames(integration) %in% seqnames(newSeqInfoData38_noextra)
      ]
      
      matchinglevels <- seqlevels(newSeqInfoData38_noextra)[seqlevels(newSeqInfoData38_noextra) %in% seqlevels(integration_noextra)]
      
      seqlengths(integration_noextra) = seqlengths(newSeqInfoData38_noextra)[seqlevels(newSeqInfoData38_noextra) %in% seqlevels(integration_noextra)]
      
      genoma_hg38 <- GRanges(
        seqnames(newSeqInfoData38_noextra),
        IRanges(
          start = 1,
          end = unname(seqlengths(newSeqInfoData38_noextra))
        )
      )
      seqinfo(genoma_hg38) <- newSeqInfoData38_noextra
      
      
      integration_plot <- Ideogram(integration_noextra, color = "red", fill = "red", alpha = 0.7,
                                       zoom.offset = 0.2, size = 1,
                                       cytobands = TRUE, aspect.ratio = 1/20
                                   )

      integration_plot
    }
    
    output$chroms <- renderPlot({
      req(input$impact)
      con <- dbConnect(RSQLite::SQLite(),
                       "/home/shared_projects/TESI/tesi_cezar_grigorean/data/test.sqlite"
      ) 
      on.exit(dbDisconnect(con), add = TRUE) 
      dbSendQuery(con, paste("SELECT *
                          FROM variants
                          WHERE impact LIKE ?"),
                  params = list(c(paste(input$impact)))
      ) %>%
        dbFetch() %>%
        dplyr::filter(grepl(paste(paste0("^", input$chr, "$"), collapse = "|"), chromosome)) %>%
        dplyr::filter(grepl(paste(paste0("^", input$chr_plot, "$"), collapse = "|"), chromosome)) %>% 
        dplyr::filter(gnomad_af <= input$gnomad |
                        is.na(gnomad_af) == TRUE) %>%
        dplyr::filter(grepl(paste0("^", input$searchgene),
                            symbol,
                            ignore.case = TRUE)) %>%
        dplyr::filter(grepl(paste(input$conseq, collapse = "|"),
                            consequence)) %>%
        dplyr::filter(grepl(paste(input$HC, collapse = "|"), lof)) %>%
        chr_information_integration()
    })
    
    output$gene <- renderPlot({
      
    })
  
    
}

shinyApp(ui, server)