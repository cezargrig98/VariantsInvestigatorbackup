#funzione con 
vcfparsing <- function(chunk, db_con){
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
  #copio il VCF modificato in database
  dbWriteTable(db_con, "variants", list_tibble, append = TRUE )
  print("uploaded")
}
