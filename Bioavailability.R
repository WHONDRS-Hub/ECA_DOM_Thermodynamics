# ==== Loading libraries =========
rm(list=ls(all=T))
# Required Libraries
library(readr)
library(dplyr)
library(ggplot2)
library(reshape2)

# Task selection
task1 <- 'Plot results'
task1a <- "OClim"

setwd('C:/Users/gara009/Downloads/MATLAB Codes 3 (1)/MATLAB Codes/')
# Input files
fticrdataFilename <- "ExperimentalData/Processed_S19S_Sediments_Water_2-2_newcode.csv"
metadataFilename <- "ExperimentalData/WHONDRS_S19S_Metadata_v3.csv"
sedrespdataFilename <- "ExperimentalData/WHONDRS_S19S_Sediment_Incubations_Respiration_Rates.csv"

# Read files
tbl_fticr <- read_csv(fticrdataFilename)
tbl_meta <- read_csv(metadataFilename)
tbl_meta = tbl_meta[-1,]
tbl_resp <- read_csv(sedrespdataFilename)

# Preprocessing
# Remove non carbon sources
tbl_fticr <- tbl_fticr %>% filter(C != 0)

# Sort tables and convert sample IDs to lower case letters
sampCol <- 39
tbl_fticr_ <- tbl_fticr[,sampCol:ncol(tbl_fticr)]
tbl_fticr_ <- tbl_fticr_[,sort(names(tbl_fticr_))]
names(tbl_fticr_) <- tolower(names(tbl_fticr_))
tbl_fticr <- cbind(tbl_fticr[,1:(sampCol-1)], tbl_fticr_)
tbl_meta$Sample_ID <- tolower(tbl_meta$Sample_ID)
tbl_resp <- tbl_resp %>% arrange(Sample_ID)
tbl_resp$Sample_ID <- tolower(tbl_resp$Sample_ID)

