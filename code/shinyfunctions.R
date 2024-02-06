library(GenomicRanges)
library(IRanges)
library(ensembldb)
library(Gviz)
library(biomaRt)
library(trackViewer)
library(plyranges)
library(DBI)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)




gene_information_integration <- function(res){
  
  res$strand[res$strand == -1] <- "-"
  res$strand[res$strand == 1] <- "+"

 
  refgenome <- readRDS("/data/refgenome_model.rds")
  
  refgenome <- refgenome %>% 
    dplyr::filter(grepl(paste0("^", res$gene),
                        gene,
                        ignore.case = TRUE)) %>% 
    dplyr::filter(is.na(exon)==FALSE & is.na(transcript)==FALSE)
  
  
  
  chr <- as.character(unique(res$chromosome))
  
  
  atrack <- AnnotationTrack(res)
  
  itrack <- IdeogramTrack(chromosome = chr, genome = "hg38")
  
  GR_res <- GRanges(seqnames = unique(res$chromosome), 
                    IRanges(start = res$start, end = res$end),
                    strand = res$strand,
                    gene_name = res$symbol, 
                    gene_id = res$gene,
                    ref = res$ref,
                    alt = res$alt,
                    conseq = res$consequence,
                    pos = res$start
                    )
  
  gtrack <- GenomeAxisTrack()
  
  variants <- AnnotationTrack(GR_res, chromosome = chr, name = "Variants",
                              strand = GR_res@strand, 
                              group = c(paste0(GR_res$pos, ".", GR_res$ref, ">", GR_res$alt)),
                              genome = "hg38",
                              fill = "red", color = "red"
                              )
  
  # refgenome_subset <- subsetByOverlaps(refgenome, GR_res)
  
  
  
  genetrack <- GeneRegionTrack(refgenome, chromosome = chr, name = as.character(res$symbol[1]),
                               trancriptAnnotation = "transcript"
                               )
  getOption("Gviz.scheme")
  scheme <- getScheme()
  scheme$GeneRegionTrack$fill <- "blue"
  scheme$GeneRegionTrack$col <- NULL
  scheme$GeneRegionTrack$transcriptAnnotation <- "transcript"
  addScheme(scheme, "myScheme")
  options(Gviz.scheme = "myScheme")
  
  
  
  plotTracks(list(itrack, gtrack, genetrack, variants), 
             background.panel = "#FFFEDB", 
             background.title = "black", groupAnnotation = "group", just.group = "left")
  }


# file_input_f <- function(){
#   directory_file <- file.choose()
#   return(dir(directory_file))
# }

#  return(dir(directory_file))
#}
