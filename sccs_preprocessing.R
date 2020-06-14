setwd("/home/jsjukara/")
library(data.table) # for using fread
library(tidyverse)
library(SCCS)
library(lubridate) # for manipulating dates
library(tictoc)


# get phenotype file
# list columns to keep from phenotype files
keepcols_end <- c("FINNGENID", "BL_AGE", "BL_YEAR", "SEX")
keepcols_endbig <- c("FINNGENID","FU_END_AGE",
                     'DEATH','DEATH_AGE','DEATH_YEAR')

# load two phenotype data files
end <- fread("zcat R5_COV_PHENO_V1.txt.gz")
end <- end[,..keepcols_end] # only keep columns of interest, done with data.table
ids <- unique(end$FINNGENID)

endbig <- fread("finngen_R5_V2_endpoint.txt")
endbig <- endbig[FINNGENID %in% ids,..keepcols_endbig]

# join small phenotype file to big phenotype file
phenotype_data <- end
phenotype_data <- left_join(phenotype_data, endbig, by="FINNGENID")

# calculate time variables
phenotype_data$YEAR_OF_BIRTH <- phenotype_data$BL_YEAR - phenotype_data$BL_AGE # year of blood sample minus age at blood sample
# calculate follow-up start age starting from 1999, since outpatient visits were added to HILMO in 1998
# but we know YEAR_OF_BIRTH only up to a year (max error is 1 year)
phenotype_data$FU_START_AGE <- 1999 - phenotype_data$YEAR_OF_BIRTH
phenotype_data$FUTIME <- phenotype_data$FU_END_AGE - phenotype_data$FU_START_AGE

# get longitudinal phenotypes
events <- fread("./finngen_data/finngen_R5_V2_endpoint_longitudinal.txt", data.table=FALSE)

events <- events[,c("FINNGENID", "EVENT_AGE", "ENDPOINT")]

# join follow-up start and end age
events <- left_join(events, phenotype_data[,c("FINNGENID", "FU_START_AGE", "FU_END_AGE")])

# filter event rows before start of follow-up
events <- events %>% filter(EVENT_AGE > FU_START_AGE)

#read longitudinal purchases
purchases <- fread("finngen_R5_v2_detailed_longitudinal.txt", data.table=FALSE)
purchases$EVENT_YEAR <- decimal_date(as.Date(purchases$APPROX_EVENT_DAY))
names(purchases)[5:8] <- c("ATC_CODE", "SAIR", "VNRO", "PLKM")

# filter rows for drug purchases after 1999 for ATC codes of interest
purchases <- purchases %>% filter(EVENT_YEAR>1999,
                                  SOURCE == "PURCH")

purchases <- purchases[,c("FINNGENID", "EVENT_AGE", "EVENT_YEAR", "ATC_CODE", "VNRO", "PLKM")]
names(purchases)[2] <- "PURCHASE_AGE"

rm(end)
rm(endbig)
rm(phenotype_data)

save.image(file = "sccs_workspace.RData", compress=FALSE)