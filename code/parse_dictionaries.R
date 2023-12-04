library(GenomicRanges)
library(IRanges)
library(ensembldb)
library(Gviz)


genome_information_integration <- function(res){
  
  
  con <- dbConnect(RSQLite::SQLite(),
                   "/home/shared_projects/TESI/tesi_cezar_grigorean/data/test.sqlite"
  )
 res <- dbSendQuery(con, paste("SELECT *
                          FROM variants
                          WHERE impact LIKE 'HIGH' ")
  ) %>%
    dbFetch() %>% 
    dplyr::filter(lof == "HC") 

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
    res$chromosome,
    IRanges(
      start = res$start,
      end = res$end
    )
  )
  
  integration_noextra <- integration[
    seqnames(integration) %in% seqnames(newSeqInfoData38_noextra)
  ]
  
  matchinglevels <- seqlevels(newSeqInfoData38_noextra)[seqlevels(newSeqInfoData38_noextra) %in% seqlevels(integration_noextra)]
  
  integration_noextra <- keepSeqlevels(integration_noextra, matchinglevels)
  
  seqlengths(integration_noextra) = seqlengths(newSeqInfoData38_noextra)[seqlevels(newSeqInfoData38_noextra) %in% seqlevels(integration_noextra)]
  
  genoma_hg38 <- GRanges(
    seqnames(newSeqInfoData38_noextra),
    IRanges(
      start = 1,
      end = unname(seqlengths(newSeqInfoData38_noextra))
    )
  )
  seqinfo(genoma_hg38) <- newSeqInfoData38_noextra
  
  ensdb <- dbConnect(RSQLite::SQLite(), "/home/shared_projects/TESI/tesi_cezar_grigorean/code/Homo_sapiens.GRCh38.110.sqlite")

genome()
genes_db <- dbSendQuery(ensdb,
                        "SELECT *
                                   FROM gene
                                   WHERE gene_id LIKE ?",
                        params = list(c(paste(res$gene, sep = ",")))
) %>% dbFetch()

genes_db <- genes_db %>% 
  mutate(gene_width = gene_seq_end - gene_seq_start + 1)
  
  genome
  atrack <- AnnotationTrack()
  plot(atrack)
  # 
  # ensdb_genes_match <- ensdb_genes_match %>%
  #   rename(gene = gene_id) %>% 
  #     desc(gene)
  # 
  # res <- full_join(res, ensdb_genes_match, by = "gene")
  # 
}

    


