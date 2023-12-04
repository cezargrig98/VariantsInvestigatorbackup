library(shiny)
library(bslib)
library(flexdashboard)
library(shinyWidgets)
library(gt)
library(tibble)
library(DT)
library(tidyverse)
# library(ggplot2)
library(ggbio)
library(plyranges)
library(GenomicRanges)

usable <- readRDS("subset_vars.rds")

ui <- page_sidebar(
  title = "VariantsInvestigator",
  sidebar = sidebar(
    title = "Filters",
    sliderInput(
      inputId = "gnomAD",
      "Maximum gnomAD allele frequency",
      min = 0,
      max = 1,
      step = 0.001,
      value = 1
    ),
    sliderInput(
      inputId = "geno_qual",
      "minimum genotype quality",
      min = 0,
      max = 99,
      step = 1,
      value = 0
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
      selected = c("HIGH", "MODERATE"),
      options = list(`actions-box` = TRUE),
      multiple = TRUE
    ),
    
    sliderInput(
      inputId = "start",
      "Section of the genome in analysis",
      min = min(usable$start),
      max = max(usable$start),
      value = c(min(usable$start), max(usable$start)),
      step = 50000
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
    
    pickerInput(
      "HC",
      "High confidence LOF variants",
      choices = c("HC"),
      selected = c(),
      multiple = TRUE
    )
    
    
  ),
  layout_columns(
    card(
      card_header("VCF_sample"),
      gt_output(outputId = "table"),
      full_screen = TRUE
    ),
    card(
      # dropdownButton(
      #
      #   tags$h3("List of Inputs"),
      #   sliderInput("mess",
      #           "Solving label overlap problem in the graphic part",
      #           min = 2,
      #           max = 10,
      #           step = 1,
      #           value = 2),
      #
      #   style = "unite", icon = icon("gear"),
      #   status = "success", width = "100px"
      #   ),
      full_screen = TRUE,
      div(style = 'height:600px; width:px; overflow: scroll',
          plotOutput("gene")),
      width = "400px"
    ),
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
    usable %>%
      dplyr::filter(SIFT_values <= input$SIFT |
                      is.na(SIFT_values) == TRUE) %>%
      dplyr::filter(PolyPhen_values >= input$PolyPhen |
                      is.na(PolyPhen_values) == TRUE) %>%
      dplyr::filter(gnomAD_AF_1 <= input$gnomAD |
                      is.na(gnomAD_AF_1) == TRUE) %>%
      dplyr::filter(GQ >= input$geno_qual)  %>%
      dplyr::filter(IMPACT_1 %in% input$impact) %>%
      dplyr::filter(start >= input$start[1] &
                      start <= input$start[2]) %>%
      dplyr::filter(grepl(input$searchgene,
                          SYMBOL_1,
                          ignore.case = TRUE)) %>%
      dplyr::filter(grepl(paste(input$conseq, collapse = "|"),
                          Consequence_1)) %>%
      dplyr::filter(grepl(paste(input$HC, collapse = "|"), LoF_1)) %>%
      dplyr::select(
        "seqnames",
        "start",
        "end",
        "strand",
        "REF",
        "ALT",
        "GQ",
        "SYMBOL_1",
        "PolyPhen_values",
        "PolyPhen_prediction",
        "SIFT_values",
        "SIFT_prediction",
        "gnomAD_AF_1",
        "Consequence_1",
        "IMPACT_1",
        "LoF_1",
        "LoF_flags_1",
        "LoF_filter_1"
      ) %>%
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
        page_size_default = 100,
        page_size_values = c(25, 50, 100, 200)
      )
  })
  
  output$gene <- renderPlot({
    usable %>%
      dplyr::filter(SIFT_values <= input$SIFT |
                      is.na(SIFT_values) == TRUE) %>%
      dplyr::filter(PolyPhen_values >= input$PolyPhen |
                      is.na(PolyPhen_values) == TRUE) %>%
      dplyr::filter(gnomAD_AF_1 <= input$gnomAD |
                      is.na(gnomAD_AF_1) == TRUE) %>%
      dplyr::filter(GQ >= input$geno_qual)  %>%
      dplyr::filter(IMPACT_1 %in% input$impact) %>%
      dplyr::filter(start >= input$start[1] &
                      start <= input$start[2]) %>%
      dplyr::filter(grepl(input$searchgene,
                          SYMBOL_1,
                          ignore.case = TRUE)) %>%
      dplyr::filter(grepl(paste(input$conseq, collapse = "|"),
                          Consequence_1)) %>%
      dplyr::filter(grepl(paste(input$HC, collapse = "|"), LoF_1)) %>%
      dplyr::select(
        "seqnames",
        "start",
        "end",
        "width",
        "strand",
        "REF",
        "ALT",
        "GQ",
        "SYMBOL_1",
        "PolyPhen_values",
        "PolyPhen_prediction",
        "SIFT_values",
        "SIFT_prediction",
        "gnomAD_AF_1",
        "Consequence_1",
        "IMPACT_1",
        "LoF_1",
        "LoF_flags_1",
        "LoF_filter_1"
      ) %>%
      dplyr::group_by(SYMBOL_1) %>%
      dplyr::arrange(start, by_group = FALSE) %>%
      dplyr::filter(SYMBOL_1 != ".") %>%
      as_granges(seqnames = seqnames, start = start, end = end, width = width, strand = strand) %>%
      autoplot(aes(color = IMPACT_1, fill = IMPACT_1)) +
        facet_wrap(~ SYMBOL_1)
      # # Ideogram(genome = "hg38")
      # ggplot(aes(
      #   x = SYMBOL_1, y = start
      #   )) +
      # geom_segment(aes(
      #   x = SYMBOL_1,
      #   xend = SYMBOL_1,
      #   y = 0,
      #   yend = start
      # )) +
      # geom_point(aes(fill=IMPACT_1),
      #   size = 3,
      #   alpha = 0.7,
      #   shape = 21,
      #   stroke = 2
      # )  +
      # theme(axis.text.x = element_text(angle = 90,
      #                                  vjust = 0.5,
      #                                  hjust = 0.5)) +
      # title("Varianti selezionate")

  })
}

shinyApp(ui, server)