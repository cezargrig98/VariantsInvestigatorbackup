library(Rsamtools)
library(dbplyr)
library(dplyr)
library(VariantAnnotation)
library(tibble)
library(DBI)
library(RSQLite)
library(ggbio)
library(plyranges)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(EnsDb.Hsapiens.v86)

ensdb_v86 <- EnsDb.Hsapiens.v86

autoplot(vcftab, type = "default", geometry = "karyogram")

makeGRangesFromDataFrame(usable, start.field = as.character(usable$start), end.field = as.character(usable$end))
?GRanges
view(seqlengths(Hsapiens))

newSeqInfoData <- Seqinfo(
  seqnames = "chr1",
  seqlengths = 248956422,
  isCircular = NA,
  genome = "Human"
)

autoplot_obj<- GRanges(usable_iranges)

autoplot(ensdb_v86, layout = "karyogram", GRangesFilter(usable_iranges), names.expr = "SYMBOL_1")

usable_iranges <- as_granges(dplyr::filter(usable, usable$IMPACT_1 == "HIGH"),
                             seqnames = seqnames, start = start, end = end, width = width, strand = STRAND_1)

usable_iranges

posinfo <- as_tibble(rowRanges(vcftab))

param <- ScanVcfParam(yieldSize=10000)

vcftab <- vcfchunking[[1]]


as_tibble(geno(vcftab)$GT)

vcfparsing <- function(vcftab){
  #raccolgo informazioni su genotipo
  GT <- as_tibble(geno(vcftab)$GT) %>%
    rename("GT" = "NG2191_NG2191")
  
  #informazioni sull'allele depth
  AD <- as_tibble(geno(vcftab)$AD) %>%
    rename("AD" = "NG2191_NG2191") %>%
    unnest_wider(col = "AD", names_sep = "_")
  #informazioni sulla depth
  DP <- as_tibble(geno(vcftab)$DP) %>%
    rename("DP" = "NG2191_NG2191")
  #informazione sulla genotype quality
  GQ <- as_tibble(geno(vcftab)$GQ) %>% 
    rename("GQ" = "NG2191_NG2191")
  
  get_alternative <- function(alt_list){
    alleles = paste0(unlist(lapply(alt_list, as.character)), collapse = ",")
    return(alleles)
  }
  get_first <- function(x){
    res = NULL
    if(is.list(x)){
      res = x[[1]]
    }
    else {
      res = x
    }
    return(x)
  }
  
  #unisco in un tibble le varie informazioni sul genotipo
  genotype <- bind_cols(GT, AD, DP, GQ)
  print("genotype")
  
  rm(GT, AD, DP, GQ)
  
  #creo un tibble con le info del VCF (precedentemente subsettate con ScanVcfParam)
  annotations_plain <- as_tibble(info(vcftab)) %>%
    select(-"CSQ") %>%
    unnest_wider(col = where(~ type_sum(.x) == "list")|where(~ type_sum(.x) == "I<list>"), names_sep = "_")
  
  print("info")
  
  #prendo l'allele ALT come colonna di stringhe per poterlo inserire su un database sql
  
  #prendo coordinate genomiche
  var_coordinates <- as_tibble(rowRanges(vcftab)) %>%
    rowwise %>%
    mutate(
      ALT <- get_alternative(ALT)
    )
    # rename("ALT"="get_alternative(ALT)")
  
  print("coordinates")
  
  
  #raccolgo informakioni sulla posizione delle varianti + creo una colonna dell'alt leggibile come char e non come DNA_string
  test_alt <- var_coordinates %>%
    filter(length(ALT) > 1) %>% 
    head(1L) %>% 
    pull(ALT)
  
  unlist(lapply(test_alt, as.character))
  
  
  
  get_alternative(test_alt)
  
  
  #unisco informazioni e posizione delle varianti
  list_tibble <- bind_cols(
    var_coordinates,
    genotype,
    annotations_plain
  ) %>% 
    mutate(
      across(where(~ type_sum(.x) == "fct"), as.character)
    # ) %>% 
    # dplyr::select(
    #   where(~ type_sum(.x) == "int") | 
    #     where(~ type_sum(.x) == "dbl") | 
    #     where(~ type_sum(.x) == "chr")
    )
  
  rm(var_coordinates,
     genotype,
     annotations_plain)
  print("unite colonne")
  #separo polyphen e sift in base a predizione di patogenicità e score puramente numerico
  list_tibble <- separate_wider_delim(list_tibble, "SIFT_1", delim = "(", names = c("SIFT_prediction","SIFT_values"), too_few = "align_start")
  list_tibble <- separate_wider_delim(list_tibble, "PolyPhen_1", delim = "(", names = c("PolyPhen_prediction","PolyPhen_values"), too_few = "align_start")
  
  
  #risistemo nome della colonna SIFT per l'app finale
  list_tibble$SIFT_values <- gsub(")", "", list_tibble$SIFT_values) %>%
    as.numeric(list_tibble$SIFT_values)
  
  #stesso lavoro per polyphen
  list_tibble$PolyPhen_values <- gsub(")", "", list_tibble$PolyPhen_values) %>%
    as.numeric(list_tibble$PolyPhen_values)  
  #sistemo colonna alt e consequence
  # usable <- usable %>%
    # rename("ALT"="ALT <- get_alternative(ALT)")
  
  # usable$Consequence_1 <- sub("(;<=\\?).*$", "", usable$Consequence_1)

  

  
  #copio il VCF modificato in database
}

usable <- list_tibble

# view(usable)
usable <- usable %>%
  dplyr::select(
      where(~ type_sum(.x) == "int") |
        where(~ type_sum(.x) == "dbl") |
        where(~ type_sum(.x) == "chr")
  ) %>% 
  rename("ALT"="ALT <- get_alternative(ALT)")
  


usable$Consequence_1 <- sub("(;<=\\?).*$", "", usable$Consequence_1)

head(table)

view(info(header(table)))

table(usable$FILTER)

as_tibble(rowRanges(table))

filt <- FilterRules(list(culpritrm = function(x) is.na(info(x)$culprit)))

half_subset <- filterVcf("whole_gen_subset.vcf", "hg19", tempfile(), filters = filt)

