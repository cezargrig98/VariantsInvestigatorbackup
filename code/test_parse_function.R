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

con <- dbConnect(RSQLite::SQLite(), "speed_test.sqlite")

vcftab <- VcfFile("/home/shared_projects/TESI/tesi_cezar_grigorean/data/whole_genome_split_clean.vcf.gz", yieldSize = 100000)
open(vcftab)

subset <- ScanVcfParam(
  info = c(
    "IMPACT",
    "SYMBOL",
    "Gene",
    "EXON",
    "INTRON",
    "Consequence",
    "SIFT",
    "PolyPhen",
    "DOMAINS",
    "gnomAD_AF",
    "CLIN_SIG",
    "PHENO",
    "LoF",
    "LoF_filter",
    "LoF_flags",
    "LoF_info"
  ),
  geno = c(
    "GT",
    "DP",
    "GQ",
    "PL"
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


# rm(vcfchunking, k)


while (nrow(vcfchunking[[k]] <- readVcf(vcftab, param = subset)) ){
  # Print progress
  cat("new lines read:", k, "\n")
  k <- k + 1
}

close(vcftab)

vcfchunking <- vcfchunking[1:k-1]

z <- length(vcfchunking)


vcfparsing <- function(chunk){
  # parallel::clusterEvalQ(cs, {library(DBI);con <- dbConnect(RSQLite::SQLite(), dbname = "variant_new.sqlit")})
  #raccolgo informazioni su genotipo
  GT <- as_tibble(geno(chunk)$GT) %>%
    rename("gt" = "NG2191_NG2191")
  
  #informazioni sulla PL
  PL <- as_tibble(geno(chunk)$PL) %>%
    rename("pl" = "NG2191_NG2191")
  PL[] <- lapply(PL, as.character)
  PL$pl <- gsub("c", "", PL$pl)
  PL$pl <- gsub(")", "", PL$pl)
  PL$pl <- sub(".", "", PL$pl)
  #informazioni sulla depth
  DP <- as_tibble(geno(chunk)$DP) %>%
    rename("dp" = "NG2191_NG2191")
  #informazione sulla genotype quality
  GQ <- as_tibble(geno(chunk)$GQ) %>% 
    rename("gq" = "NG2191_NG2191")
  
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
  genotype <- bind_cols(gt, pl, dp, gq)
  print("genotype")
  
  rm(GT, PL, DP, GQ)
  
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
    dplyr::filter(length(ALT) > 1) %>% 
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
  # dbWriteTable(con, "variants", list_tibble, append = TRUE)
  # print("uploaded")
  
  #copio il VCF modificato in database
}




registerDoParallel(cores = 15)


foreach(i=1:z) %dopar% {
  vcfparsing(vcfchunking[[i]])
} 

ready_to_go <- bind_cols(vcfchunking)

rm(vcfchunking)

copy_to(con, ready_to_go, "variants",
        temporary = FALSE, 
        indexes = list(
          c("seqnames", "start", "end", "ref", "ALT"))
)

rm(ready_to_go)
