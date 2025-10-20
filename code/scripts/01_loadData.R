# Introduction ------------------------------------------------------------

### -- this script serves the purpose of loading our bulk-RNA seq data and exploring aspects of it. 
### -- this data was obtained from Laura of Jen Jen Yeh's lab, and stored in:
  ## -- https://adminliveunc-my.sharepoint.com/:f:/r/personal/xlpeng_ad_unc_edu/Documents/Collaborations/MDACC_evRNA/PDAC?csf=1&web=1&e=rZ4GNf

rm(list = ls())
options(timeout = 1e6)

# Load libraries ----------------------------------------------------------

library(tidyverse)
library(glue)
library(magrittr)
library(DESeq2)
library(SummarizedExperiment)
library(tools)
library(glue)
library(biomaRt)


# Loading data ------------------------------------------------------------

bulkRDS_FF <- readRDS("data/02_psrc_maitra_BulkRNA_refSeq_withsubtypeCal.rds")$FF
bulkRDS_FFPE <- readRDS("data/02_psrc_maitra_BulkRNA_refSeq_withsubtypeCal.rds")$FFPE

  ## -- we need to use both the FF and FFPE data

## -- convert to a summarized experiment

se <- SummarizedExperiment(bulkRDS_FF[["counts"]])
names(assays(se)) <- "counts"

# Exploratory data analysis (EDA) -----------------------------------------

  ## -- extracting the count matrix

counts_FF <- bulkRDS_FF[["counts"]]
counts_FFPE <- bulkRDS_FFPE[["counts"]]


  ## -- extracting the evRNA seq data

### -- first extract the gene count data

ev_samps_FF <- bulkRDS_FF[["sampInfo"]] %>% filter(`Plasma/Tissue` == "Plasma") %>% rownames()
ev_samps_FFPE <- bulkRDS_FFPE[["sampInfo"]] %>% filter(`Plasma/Tissue` == "Plasma") %>% rownames()

counts_ev_FF <- counts_FF %>% select(all_of(ev_samps_FF))
counts_ev_FFPE <- counts_FFPE %>% select(all_of(ev_samps_FFPE))

  ## -- extracting only bulk tumor gene exp data

counts_bulk_FF <- counts_FF %>% select(!all_of(ev_samps_FF))
counts_bulk_FFPE <- counts_FFPE %>% select(!all_of(ev_samps_FFPE))


#### -- creating comprehensive evRNA & bulk datasets

counts_bulk <-
  counts_bulk_FF %>% 
    bind_cols(counts_bulk_FFPE)

counts_ev <-
  counts_ev_FF %>%
    bind_cols(counts_ev_FFPE)

data <- list(bulk = counts_bulk, ev = counts_ev)
saveRDS(data, "data/data_summarized.RDS")


# Public Datasets ---------------------------------------------------------


## (1) TCGA-PAAD: tried to get it but was not successful----




## (2) CPTAC ----

cptac <- read.table("data/public/CPTAC/mRNA_RSEM_UQ_log2_Tumor.cct",
                  sep = "\t", header = TRUE, row.names = 1)


## (3) Moffitt ----

moffitt <- read.table("data/public/Moffitt/GSE71729_series_matrix.txt",
                    sep = "\t", header = TRUE, comment.char = "!", row.names = 1)

  ## -- need to subset to only primary PDAC or stroma LCM samples

lines <- readLines("data/public/Moffitt/GSE71729_series_matrix.txt")

sample_ids <- grep("^!Series_sample_id", lines, value = TRUE)
sample_ids <- sub('^!Series_sample_id\\t\\"', "", sample_ids)
sample_ids <- sub('^\t\"', "", sample_ids); sample_ids <- sub('^\"', "", sample_ids)

sample_source <- grep("^!Sample_source_name_ch2", lines, value = TRUE)
sample_source <- sub('^!Sample_source_name_ch2\\t\\"', "", sample_source)

sample_ids <- unlist(strsplit(sample_ids, " "))
sample_ids <- sample_ids[which(sample_ids %in% grep("^GSM", sample_ids, value = TRUE))]

sample_source <- unlist(strsplit(sample_source, '"\\t\\"'))

ids_primary <- sample_ids[which(sample_source %in% "Pancreas_Primary")]

  ## -- do the subsetting

moffitt %<>%
  dplyr::select(all_of(ids_primary))



## (4) DIJK ----

dijk <- read.table("data/public/DIJK/E-MTAB-6830.sdrf.txt",
                   sep = "\t", header = TRUE)

fun_download <- function(x) {
  
  file_name = file_path_sans_ext(file_path_sans_ext(basename(x)))
  
  download.file(x, 
                destfile = 
                  glue("data/public/DIJK/ind_files/{file_name}.gz"), 
                       method = "auto")
}

sapply(X = dijk$Comment.FASTQ_URI., fun_download)

  ## -- opening a single file and checking

df <- read.delim("data/public/DIJK/ind_files/ERR2603466.gz")


## (5) Linehan ----

linehan <- read.table(file = "data/public/LINEHAN/GSE131050_PurIST_Linehan_seq.tsv",
                      sep = "\t", header = TRUE, row.names = 1)

rownames(linehan) <- sub("^[^;]*;", "", rownames(linehan))

linehan %<>%
  tibble::rownames_to_column(var = "gene") %>%
  arrange(gene) %>%
  tibble::column_to_rownames(var = "gene")


## (6) PACA AU-seq ----



## (7) PACA AU-array ----



## (8) Puleo array ----



## (9) Hayashi (RNA) ----



## (10) Sears ---- this is a more recent dataset

sears <- read.table("data/public/Sears/GSE205154_Gene_Level_TPM_Estimates.txt",
                    sep = "\t", header = TRUE)






