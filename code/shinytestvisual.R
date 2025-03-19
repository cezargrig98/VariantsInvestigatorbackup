library(GenomicRanges)
library(IRanges)
library(Gviz)
library(trackViewer)
library(DBI)
library(shiny)
library(shiny.gosling)
library(RSQLite)

function(res){
  con <- dbConnect(RSQLite::SQLite(), "data/4035-23N.haplotypecaller_hardfiltered_PASS_ontarget_vepped_split_cleansql.sqlite")
  
  res <- dbGetQuery(con, 
                    "SELECT * 
                    FROM variants
                    WHERE impact LIKE 'HIGH'
                    ") |>
    dplyr::filter(symbol == "STXBP2")
  
  res$strand[res$strand == -1] <- "-"
  res$strand[res$strand == 1] <- "+"
  
  
  refgenome <- readRDS("data/refgenome_model.rds")
  
  refgenome <- refgenome |> 
    dplyr::filter(grepl(paste0("^", 	
                               res$gene),
                        gene,
                        ignore.case = TRUE)) |> 
    dplyr::filter(is.na(exon)==FALSE & is.na(transcript)==FALSE) 
  # |> 
  #   dplyr::filter(start >= res$start && end <= res$end)
  
  GR_res <- GRanges(seqnames = unique(res$chromosome), 
                    IRanges(start = res$start, end = res$end, names = res$id)
                    
  )
  
  features <- GRanges(seqnames = unique(refgenome$chromosome), 
                      IRanges(start = refgenome$start, end = refgenome$end, names = refgenome$transcript)
  )
  features <- subsetByOverlaps(features, GR_res)
  features$fill <- c("#FF8833", "#51C6E6")
  
  features.mul$height[1] <- list(unit(1/8, "inches"),
                                 unit(0.5, "lines"),
                                 unit(.2, "char"))
  features.mul$height[2] <- list(unit(2/8, "inches"),
                                 unit(1, "lines"),
                                 unit(.4, "char"))
  
  features.mul$featureLayerID <- features.mul$id
  
  names(features.mul) <- features.mul$featureLayerID
  
  lolliplot(GR_res, features.mul)
}