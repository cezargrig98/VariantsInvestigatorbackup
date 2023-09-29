
#R code stuff-----------
#carico i pacchetti necessari al funzionamento dell'app
library(VariantAnnotation)
library(tibble)
library(Rsamtools)
library(dbplyr)
library(RSQLite)
library(DBI)
library(tidyr)
library(flexdashboard)
library(shiny)
library(purrr)
library(parallel)


indexTabix("/home/shared_projects/TESI/tesi_cezar_grigorean/data/whole_genome_split.vcf.gz.gz", "vcf")
#creo database sql per il VCF
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "new_variants.sqlit")
## parte di apertura del VCF

vcfchunking <- readVcf(TabixFile("/home/shared_projects/TESI/tesi_cezar_grigorean/data/whole_genome_split.vcf.gz.gz", yieldSize = 1000))
## parte in cui unisco tutto il vcf tramite ciclo while
vcf <- list()
k <- 1
while (nrow( vcf[[k]] <- readVcf(vcfchunking))){
  print(k)
  
  #raccolgo informazioni sulle varianti
  GT <- as.tibble(geno(vcf[[z]])$GT) %>%
    rename("GT" = "NG2191_NG2191")
  
  
  AD <- as.tibble(geno(vcf[[z]])$AD) %>%
    rename("AD" = "NG2191_NG2191") %>%
    unnest_wider(col = "AD", names_sep = "_")
  
  DP <- as_tibble(geno(vcf[[z]])$DP) %>%
    rename("DP" = "NG2191_NG2191") %>%
    unnest_wider(col="DP", names_sep = "_")
  
  GQ <- as_tibble(geno(vcf[[z]])$GQ) %>% 
    rename("GQ" = "NG2191_NG2191")
  
  
  MIN_DP <- as_tibble(geno(vcf[[z]])$MIN_DP) %>% 
    rename("MIN_DP" = "NG2191_NG2191")
  
  
  PGT <- as_tibble(geno(vcf[[z]])$PGT) %>% 
    rename("PGT" = "NG2191_NG2191")
  
  
  PID <- as_tibble(geno(vcf[[z]])$PID) %>% 
    rename("PID" = "NG2191_NG2191")
  
  
  PL <- as_tibble(geno(vcf[[z]])$PL) %>%
    rename("PL" = "NG2191_NG2191") %>%
    unnest_wider(col = "PL", names_sep = "_")
  
  PS <- as_tibble(geno(vcf[[z]])$PS) %>% 
    rename("PS" = "NG2191_NG2191")
  
  RGQ <- as_tibble(geno(vcf[[z]])$RGQ) %>% 
    rename("RGQ" = "NG2191_NG2191")
  
  
  
  genotype <- bind_cols(GT, AD, DP, GQ, MIN_DP, PGT, PID, PL, PS, RGQ)
  
  
  
  annotations_plain = as_tibble(info(vcf[[z]]))%>%
    unnest_wider(col = where(~ type_sum(.x) == "list")|where(~ type_sum(.x) == "I<list>"), names_sep = "_")
  
  
  
  
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
  var_coordinates = as_tibble(rowRanges(vcf[[z]])) %>%
    rowwise %>%
    mutate(
      ALT = get_alternative(ALT)
    ) 
  
  
  
  
  #rimuovo dalle info le colonne superflue, rendo digeribili le info al tibble
  annotations_plain = as_tibble(info(vcf[[z]]))%>%
    unnest_wider(col = where(~ type_sum(.x) == "list")|where(~ type_sum(.x) == "I<list>"), names_sep = "_")
  
  #raccolgo informazioni sulla posizione delle varianti
  test_alt = var_coordinates %>%
    filter(length(ALT) > 1) %>% 
    head(1L) %>% 
    pull(ALT)
  
  unlist(lapply(test_alt, as.character))
  
  
  
  get_alternative(test_alt)
  
  
  #unisco informazioni e posizione delle varianti
  
  vcf[[z]] = bind_cols(
    var_coordinates,
    genotype,
    annotations_plain
  ) %>% 
    dplyr::select(-("END")) %>% 
    mutate(
      across(where(~ type_sum(.x) == "fct"), as.character)
    ) %>% 
    dplyr::select(
      where(~ type_sum(.x) == "int") | 
        where(~ type_sum(.x) == "dbl") | 
        where(~ type_sum(.x) == "chr")
    )
  
  #separo polyphen e sift in base a predizione di patogenicità e score puramente numerico
  vcf[[z]] <- separate_wider_delim(vcf[[z]], "SIFT_1", delim = "(", names = c("SIFT_prediction","SIFT_values"), too_few = "align_start")
  vcf[[z]] <- separate_wider_delim(vcf[[z]], "PolyPhen_1", delim = "(", names = c("PolyPhen_prediction","PolyPhen_values"), too_few = "align_start")
  
  
  vcf[[z]]$SIFT_values <- gsub(")", "", vcf[[z]]$SIFT_values) %>%
    as.numeric(vcf[[z]]$SIFT_values)
  vcf[[z]]$PolyPhen_values <- gsub(")", "", vcf[[z]]$PolyPhen_values) %>%
    as.numeric(vcf[[z]]$PolyPhen_values)
  
  rm(GT,AD,DP,GQ,MIN_DP,PGT,PID,PL,PS,RGQ)
  
  k <- k + 1
  
}

close(vcfchunking)

#elimino vcf vuoto
vcf <- vcf[1:k-1]

#modifico in tibble tutti i chunk del VCF
z <- 1



while (z < k) {
#raccolgo informazioni sulle varianti
GT <- as.tibble(geno(vcf[[z]])$GT) %>%
  rename("GT" = "NG2191_NG2191")
  

AD <- as.tibble(geno(vcf[[z]])$AD) %>%
  rename("AD" = "NG2191_NG2191") %>%
  unnest_wider(col = "AD", names_sep = "_")

DP <- as_tibble(geno(vcf[[z]])$DP) %>%
  rename("DP" = "NG2191_NG2191") %>%
  unnest_wider(col="DP", names_sep = "_")

GQ <- as_tibble(geno(vcf[[z]])$GQ) %>% 
  rename("GQ" = "NG2191_NG2191")


MIN_DP <- as_tibble(geno(vcf[[z]])$MIN_DP) %>% 
  rename("MIN_DP" = "NG2191_NG2191")


PGT <- as_tibble(geno(vcf[[z]])$PGT) %>% 
  rename("PGT" = "NG2191_NG2191")


PID <- as_tibble(geno(vcf[[z]])$PID) %>% 
  rename("PID" = "NG2191_NG2191")


PL <- as_tibble(geno(vcf[[z]])$PL) %>%
  rename("PL" = "NG2191_NG2191") %>%
  unnest_wider(col = "PL", names_sep = "_")

PS <- as_tibble(geno(vcf[[z]])$PS) %>% 
  rename("PS" = "NG2191_NG2191")

RGQ <- as_tibble(geno(vcf[[z]])$RGQ) %>% 
  rename("RGQ" = "NG2191_NG2191")



genotype <- bind_cols(GT, AD, DP, GQ, MIN_DP, PGT, PID, PL, PS, RGQ)


  
annotations_plain = as_tibble(info(vcf[[z]]))%>%
  unnest_wider(col = where(~ type_sum(.x) == "list")|where(~ type_sum(.x) == "I<list>"), names_sep = "_")


  
  
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
  var_coordinates = as_tibble(rowRanges(vcf[[z]])) %>%
    rowwise %>%
    mutate(
      ALT = get_alternative(ALT)
      ) 

 
  
  
  #rimuovo dalle info le colonne superflue, rendo digeribili le info al tibble
  annotations_plain = as_tibble(info(vcf[[z]]))%>%
    unnest_wider(col = where(~ type_sum(.x) == "list")|where(~ type_sum(.x) == "I<list>"), names_sep = "_")
  
#raccolgo informazioni sulla posizione delle varianti
test_alt = var_coordinates %>%
    filter(length(ALT) > 1) %>% 
    head(1L) %>% 
    pull(ALT)

unlist(lapply(test_alt, as.character))



get_alternative(test_alt)


#unisco informazioni e posizione delle varianti

vcf[[z]] = bind_cols(
    var_coordinates,
    genotype,
    annotations_plain
    ) %>% 
    dplyr::select(-("END")) %>% 
    mutate(
      across(where(~ type_sum(.x) == "fct"), as.character)
    ) %>% 
    dplyr::select(
      where(~ type_sum(.x) == "int") | 
        where(~ type_sum(.x) == "dbl") | 
        where(~ type_sum(.x) == "chr")
    )

#separo polyphen e sift in base a predizione di patogenicità e score puramente numerico
vcf[[z]] <- separate_wider_delim(vcf[[z]], "SIFT_1", delim = "(", names = c("SIFT_prediction","SIFT_values"), too_few = "align_start")
vcf[[z]] <- separate_wider_delim(vcf[[z]], "PolyPhen_1", delim = "(", names = c("PolyPhen_prediction","PolyPhen_values"), too_few = "align_start")


vcf[[z]]$SIFT_values <- gsub(")", "", vcf[[z]]$SIFT_values) %>%
  as.numeric(vcf[[z]]$SIFT_values)
vcf[[z]]$PolyPhen_values <- gsub(")", "", vcf[[z]]$PolyPhen_values) %>%
  as.numeric(vcf[[z]]$PolyPhen_values)

rm(GT,AD,DP,GQ,MIN_DP,PGT,PID,PL,PS,RGQ)

  z <- z+1
}

#unisco pezzi del vcf
vcf <- bind_rows(vcf)%>%
  select(-("LOF_1"))

#carico vcf in sqlite
copy_to(con, vcf, "variants",
        temporary = FALSE, 
        indexes = list(
          c("seqnames", "start", "end", "ref", "ALT"))
        )
#pulizia
rm(z,k, annotations_plain,vcf,var_coordinates,get_alternative,get_first,vcf,chunksize,test_alt,vcfchunking,genotype)

variants_db <- tbl(con, "variants")

view(variants_db)
#---------------

#funzioni filtro-------------
#frequenza gnomAD
freq_gnomAD <- function(gnomin, gnomax){variants_db %>%
    dplyr::filter(between(gnomAD_AF_1, gnomAD[1], gnomAD[2]))
}
#frequenza massima tra gnomAD, 1000genome e ESP
freq_max <- function(min_maf, max_maf){variants_db %>%
  dplyr::filter(between(MAX_AF, min_maf, max_maf))
  }
#frequenza 1000genome
freq_1000genome <- function(min1000, max1000){variants_db %>%
  dplyr::filter(between(AF, min1000, max1000))
  }
#genotype quality
geno_qual <- function(gq){variants_db %>%
  dplyr::filter(GenoQual >= gq)
  }
#read depth
read_depth <-function(read_dp){variants_db %>%
  dplyr::filter(readDP >= read_dp)
  }
#gene name
gene_name <- function(name){variants_db %>%
  dplyr::filter(SYMBOL_1 == name)
}
#dove cade la mutazione (conseguenza) -> da
var_conseq <- function(conseq){variants_db %>%
  dplyr::filter(Consequence %like% conseq)
  }
#esone/introne della mutazione
exon <- function(ex){variants_db %>%
  dplyr::filter(EXON == ex)
  }
intron <- function(intr){variants_db %>%
  dplyr::filter(INTRON == intr)
}
#SIFT e PolyPhen
sift <- function(minsft, maxsft){variants_db %>%
    dplyr::filter(between(SIFT_values, minsft, maxsft))
}
polyphen <- function(minpol, maxpol){variants_db %>%
    dplyr::filter(between(PolyPhen_values, minpol, maxpol))
}
#flag per lof annotate
flag_lof <- function(){variants_db %>%
    filter(LoF_1 == "HC")}

polyphen(0,1)

view(variants_db %>% select(LoF_1))
#----------------------------

#------------------------------------------