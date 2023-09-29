##chunking vcf e subset delle informazioni più importanti
library(VariantAnnotation)
library(tibble)
library(Rsamtools)
library(dbplyr)
library(dplyr)
library(RSQLite)
library(tidyr)
library(flexdashboard)
library(shiny)
library(purrr)
library(parallel)
library(doParallel)
library(DBI)
library(foreach)

#creo file index per velocizzare la lettura del vcf


## parte di apertura del VCF

indexTabix("/home/shared_projects/TESI/tesi_cezar_grigorean/data/whole_genome_split.vcf.gz", "vcf")

vcftab <- VcfFile("/home/shared_projects/TESI/tesi_cezar_grigorean/data/whole_genome_split.vcf.gz", yieldSize = 100000)
open(vcftab)

subset <- ScanVcfParam(
  info = c(
    "IMPACT",
    "SYMBOL",
    "Consequence",
    "SIFT",
    "PolyPhen",
    "gnomAD_AF",
    "LoF",
    "LoF_filter",
    "LoF_flags",
    "LoF_info"
  ),
  geno = c(
    "AD",
    "GT",
    "DP",
    "GQ"
  )
)


#fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#param <- ScanVcfParam(fixed="ALT", geno=c("GT", "GL"), info=c("LDAF"))
#tab <- VcfFile(fl, yieldSize=4000)
#open(tab)
#while (nrow(vcf_yield <- readVcf(tab, "hg19", param=param)))
#  cat("vcf dim:", dim(vcf_yield), "\n")
#close(tab)
vcfchunking <- list()


k <- 1





while (nrow(vcfchunking[[k]] <- readVcf(vcftab), param = subset) ){
  # Print progress
  cat("new lines read:", k, "\n")
  k <- k + 1
}

vcfchunking <- vcfchunking[1:k-1]

z <- length(vcfchunking)


vcfparsing <- function(chunk){
  parallel::clusterEvalQ(cs, {library(DBI);con <- dbConnect(RSQLite::SQLite(), dbname = "variant_new.sqlit")})
  #raccolgo informazioni su genotipo
  GT <- as_tibble(geno(chunk)$GT) %>%
    rename("GT" = "NG2191_NG2191")
  
  #informazioni sull'allele depth
  AD <- as_tibble(geno(chunk)$AD) %>%
    rename("AD" = "NG2191_NG2191") %>%
    unnest_wider(col = "AD", names_sep = "_")
  #informazioni sulla depth
  DP <- as_tibble(geno(chunk)$DP) %>%
    rename("DP" = "NG2191_NG2191")
  #informazione sulla genotype quality
  GQ <- as_tibble(geno(chunk)$GQ) %>% 
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
  annotations_plain <- as_tibble(info(chunk)) %>%
    unnest_wider(col = where(~ type_sum(.x) == "list")|where(~ type_sum(.x) == "I<list>"), names_sep = "_")
  
  print("info")
  
  #prendo l'allele ALT come colonna di stringhe per poterlo inserire su un database sql
  
  #prendo coordinate genomiche
  var_coordinates <- as_tibble(rowRanges(chunk)) %>%
    rowwise %>%
    mutate(
      ALT <- get_alternative(ALT)
    ) 
  
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
    ) %>% 
    dplyr::select(
      where(~ type_sum(.x) == "int") | 
        where(~ type_sum(.x) == "dbl") | 
        where(~ type_sum(.x) == "chr")
    )
  
  rm(var_coordinates,
     genotype,
     annotations_plain)
  print("unite colonne")
  #separo polyphen e sift in base a predizione di patogenicità e score puramente numerico
  list_tibble <- separate_wider_delim(list_tibble, "SIFT_1", delim = "(", names = c("SIFT_prediction","SIFT_values"), too_few = "align_start")
  list_tibble <- separate_wider_delim(list_tibble, "PolyPhen_1", delim = "(", names = c("PolyPhen_prediction","PolyPhen_values"), too_few = "align_start")
  print("separati sift e polyphen in prediction+valori")
  
  
  #risistemo nome della colonna SIFT per l'app finale
  list_tibble$SIFT_values <- gsub(")", "", list_tibble$SIFT_values) %>%
    as.numeric(list_tibble$SIFT_values)
  
  #stesso lavoro per polyphen
  list_tibble$PolyPhen_values <- gsub(")", "", list_tibble$PolyPhen_values) %>%
    as.numeric(list_tibble$PolyPhen_values)
  print("almost uploaded")
  dbWriteTable(con, "variants", list_tibble, append = TRUE)
  print("uploaded")
  
  #copio il VCF modificato in database
}



registerDoParallel(cores = 10)

list_tibble <- list()

foreach(i=1:z) %dopar% {
  vcfparsing(vcfchunking[[i]])
} 




