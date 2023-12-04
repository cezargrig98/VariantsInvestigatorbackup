##chunking vcf e subset delle informazioni più importanti
library(VariantAnnotation)
library(tibble)
library(Rsamtools)
library(dbplyr)
library(dplyr)
library(RSQLite)
library(DBI)
library(tidyr)
library(flexdashboard)
library(shiny)
library(purrr)
library(parallel)
library(doParallel)
library(DBI)





#creo database sql per il VCF
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = "speed_test")
## parte di apertura del VCF

#dbRemoveTable(con, "variants")

subset <- ScanVcfParam(info = c("IMPACT",
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
                       geno = c("AD",
                                "GT",
                                "DP",
                                "GQ")
                       )

vcftab <- VcfFile("/home/shared_projects/TESI/tesi_cezar_grigorean/data/whole_genome_split.vcf.gz", yieldSize = 100000)
open(vcftab)


#fl <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
#param <- ScanVcfParam(fixed="ALT", geno=c("GT", "GL"), info=c("LDAF"))
#tab <- VcfFile(fl, yieldSize=4000)
#open(tab)
#while (nrow(vcf_yield <- readVcf(tab, "hg19", param=param)))
#  cat("vcf dim:", dim(vcf_yield), "\n")
#close(tab)
vcfchunked<- list()


k <- 1



while (nrow( vcfchunked[[k]] <- readVcf(vcftab)) ){
  # Print progress
  cat("new lines read:", k, "\n")
  k <- k + 1
}

vcfchunking <- vcfchunking[1:k-1]

z <- length(vcfchunking)


registerDoParallel(15)

foreach(i=1:z) %dopar% {
  vcfparsing(vcfchunked[[i]], con)
}  

copy_to(con, vcf, "variants",
        temporary = FALSE, 
        indexes = list(
          c("seqnames", "start", "end", "ref", "ALT"))
)

variants_db <- tbl(con, "variants")

variants_db
