library(GenomicRanges)
library(IRanges)
library(Gviz)
library(trackViewer)
library(plyranges)
library(DBI)
# library(geneviewer)



#visualization function
gene_information_integration <- function(res){
  
  res$strand[res$strand == -1] <- "-"
  res$strand[res$strand == 1] <- "+"
  
  
  refgenome <- readRDS("data/refgenome_model.rds")
  
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
                    pos = res$start,
                    id = res$id
  )
  
  gtrack <- GenomeAxisTrack()
  
  variants <- AnnotationTrack(GR_res, chromosome = chr, name = "Variants",
                              strand = GR_res@strand, 
                              group = GR_res$id,
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
  
  
  
  plotTracks(list(genetrack, variants), 
             background.panel = "#FFFEDB", 
             background.title = "black", groupAnnotation = "group", just.group = "left")
}


