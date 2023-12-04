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
library(foreach)

x <- detectCores()

con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "speed_test")


size = 10000


readVcf(TabixFile(fl, yieldSize))

scanVcfHeader(myVcfFile)

myVcfFile <- VcfFile("/home/shared_projects/TESI/tesi_cezar_grigorean/data/whole_genome_split.vcf.gz", yieldSize = size)

open(myVcfFile)

chunksvcf <- list()
k <- 1

while (nrow( chunksvcf[[k]] <- readVcf(myVcfFile))){
  print(k)
  k <- k+1
}

info(chunksvcf[[1]])
info(chunksvcf[[2]])

i = 2




