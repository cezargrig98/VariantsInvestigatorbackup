

## potenziali colonne da eliminare
-("Allele"),
-("PUBMED"),
-("AFR_AF"), -("AMR_AF"), -("EAS_AF"), -("SAS_AF"), -("AA_AF"), -("EA_AF"),
-("gnomAD_AFR_AF"), -("gnomAD_AMR_AF"), -("gnomAD_ASJ_AF"), -("gnomAD_EAS_AF"), -("gnomAD_FIN_AF"), -("gnomAD_OTH_AF"), -("gnomAD_SAS_AF"),
-("DISTANCE"),
-("MANE_SELECT"), -("MANE_PLUS_CLINICAL"),
-("TREMBL"), -("SWISSPROT"), -("UNIPARC"), -("UNIPROT_ISOFORM"),
-("miRNA"),
-("MOTIF_NAME"),
-("MOTIF_POS"),
-("HIGH_INF_POS"),
-("MOTIF_SCORE_CHANGE"),
-("TRANSCRIPTION_FACTORS"),
-("TSL"),
-("APPRIS"),
-("EUR_AF"),
-("MAX_AF_POPS"),
-("gnomAD_NFE_AF"),
-("MAX_AF"),
-("LOF"),
-("InbreedingCoeff"),
-("MLEAC"),
-("MLEAF"),
-("MQ"),
-("MQRankSum"),
-("RAW_MQandDP"),
-("ReadPosRankSum"),
-("SOR"),
-("ANN"),
-("NMD"),
-("AC"),
-("BaseQRankSum"),
-("FS")
)



##funzioni per riadattare il vcf ai miei scopi
library(tidyverse)
library(VariantAnnotation)
​
vcf <- readVcf("dataset_bsa_vepped_split_clean.vcf.gz", genome = "hg38")
​
annotations <- as_tibble(info(vcf))
​
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
​
​
​
annotations_plain = annotations %>% 
  unnest_wider(col = where(~ type_sum(.x) == "list"), names_sep = "_")
​
​
var_coordinates = as_tibble(rowRanges(vcf))
​
test_alt = var_coordinates %>%
  rowwise() %>%
  filter(length(ALT) > 1) %>% 
  head(1L) %>% 
  pull(ALT)
​
unlist(lapply(test_alt, as.character))
​
​
​
get_alternative <- function(alt_list){
  alleles = paste0(unlist(lapply(alt_list, as.character)), collapse = ",")
  return(alleles)
}
​
get_alternative(test_alt)
​
​
var_coordinates = as_tibble(rowRanges(vcf)) %>% 
  rowwise() %>% 
  mutate(
    ALT = get_alternative(ALT)
  )
​
​

variants = bind_cols(
  var_coordinates,
  annotations_plain
) %>% 
  dplyr::select(-c(END, LOF_1)) %>% 
  mutate(
    across(where(~ type_sum(.x) == "fct"), as.character)
  ) %>% 
  dplyr::select(
    where(~ type_sum(.x) == "int") | 
      where(~ type_sum(.x) == "dbl") | 
      where(~ type_sum(.x) == "chr")
  )
​
library(dplyr)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "test_variants_db.sqlit")
​
​
copy_to(con, variants, "revariants",
        temporary = FALSE, 
        indexes = list(
          c("seqnames", "start", "end", "ref", "ALT")
        )
)
​
​
rm(vcf, annotations, annotations_plain, var_coordinates, variants)
​
variants_db <- tbl(con, "variants")
​
variants_db %>% filter(
  QUAL > 200
)
​

so_consequences - readRDS(url("<https://github.com/lescai-teaching/datasets_class/blob/master/reference/annotations/ranked_ensembl_consequences.RData?raw=true", "rb"))

get_from_ann <- function(annotation, element){
  data <- unlist(str_split(annotation, "\\|"))
  return(data[element])
}

get_rank <- function(consequence){
  rank = so_consequences$rank[which(so_consequences$SO_term %in% consequence)]
  return(rank)
}

get_most_severe_index <- function(annotations_list){
  consequences = unlist(lapply(annotations_list, get_from_ann, element = 2))
  ranks = unlist(lapply(consequences, get_rank))
  most_severe = which(min(ranks) %in% ranks)
  return(most_severe[1])
}

get_most_severe_consequence <- function(annotations_list){
  consequence <- tryCatch(
    {
      index <- get_most_severe_index(annotations_list)
      get_from_ann(annotations_list[[index]], element = 2)
    },
    error=function(cond){
      return(NA)
    }
  )
  return(consequence)
}

get_most_severe_gene <- function(annotations_list){
  gene <- tryCatch(
    {
      index <- get_most_severe_index(annotations_list)
      get_from_ann(annotations_list[[index]], element = 4)
    },
    error=function(cond){
      return(NA)
    }
  )
  return(gene)
}

get_most_severe_impact <- function(annotations_list){
  impact <- tryCatch(
    {
      index <- get_most_severe_index(annotations_list)
      get_from_ann(annotations_list[[index]], element = 3)
    },
    error=function(cond){
      return(NA)
    }
  )
  return(impact)
}






#geno variants

DP <- as.numeric(geno(fullvcf[[z]])$DP)
GT <- as.character(geno(fullvcf[[z]])$GT)
GQ <- as.numeric(geno(fullvcf[[z]])$GQ)
aldepth_geno <- as_tibble(geno(fullvcf[[z]])$AD)

