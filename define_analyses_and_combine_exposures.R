setwd("/home/jsjukara/")
library(data.table) # for using fread
library(tidyverse)
library(SCCS)
library(lubridate) # for manipulating dates
library(tictoc)

# load workspace which contains "events" and "purchases" data frames from preprocessing step
load("sccs_workspace.RData")

## Define endpoints and ATC-codes for different rounds of analyses

# function to go from a vector into string, e.g. c("A", "B") becomes string "A|B"
vec2sep <- function(x, sep="|") {
    return(paste(x, collapse=sep))
}

# function that reverses the above
sep2vec <- function(x, sep="|") {
    return(unlist(strsplit(x, split="\\|")))
}    
    
bleed_endpoints <- c("D3_HAEMORRHAGCIRGUANTICO", #Haemorrhagic disorder due to circulating anticoagulants
                    "H7_CONJUHAEMOR", #conjunctival haemorrhage
                    "H7_RETINAHAEMORR", #Retinal haemorrhage
                    "H7_VITRHAEMORR", #Vitreous haemorrhage
                    "I9_INTRACRA", #Nontraumatic intracranial haemmorrhage: I9_SAH|I9_ICH
                     "I9_OTHINTRACRA", #Other intracranial haemorrhages
                     "ST19_EPIDU_HAEMORRHAGE", #Epidural haemorrhage
                    "ST19_TRAUMAT_SUBDU_HAEMORRHAGE", #Traumatic subdural haemorrhage
                    "ST19_TRAUMAT_SUBAR_HAEMORRHAGE", #Traumatic subarachnoid haemorrhage
                    "R18_HAEMORRHAGE_RESPI_PASSA") #Haemorrhage from respiratory passages

myopathy_endpoints <- c("G6_DRUGMYOP", #only 43 cases
                        "M13_MYALGIA")

cvd_endpoints <- c("I9_MI", # Myocardial infarction
               "I9_STR_EXH") # Stroke, excluding SAH

bzd_endpoints <- c("ST19_FRACT_FEMUR", # Hip fractures
                   "ST19_FRACT_WRIST_HAND_LEVEL", # Wrist/hand fractures
                   "FALLS", # Falls/tendency to fall
                   "VWXY20_UNSPE_FALL",
                   "VWXY20_TRANSPO_ACCIDENTS" # Transport accidents
                  )

opioid_endpoints <- c("K11_CONSTIPATION")

agranulocytosis_endpoints <- c("D3_AGRANULOCYTOSIS")

osteoporosis_endpoints <- c("M13_OSTEOPOROSIS",
                           "ST19_FRACT_RIBS_STERNUM_THORACIC_SPINE",
                           "ST19_FRACT_LUMBAR_SPINE_PELVIS",
                           "ST19_FRACT_SHOUL_UPPER_ARM",
                           "ST19_FRACT_FOREA",
                           "ST19_FRACT_WRIST_HAND_LEVEL",
                           "ST19_FRACT_LOWER_LEG_INCLU_ANKLE",
                           "ST19_FRACT_SPINE_LEVEL_UNSPE")

clot_endpoints <- c("I9_DVTANDPULM") # DVT and PE

atc_endpoint_map <- data.frame(drug_group = c('ssri',
                                     'statin',
                                     "cox-2 inhibitor",
                                     "benzodiazepine",
                                     "opioid",
                                     "antiepileptics/clozapine",
                                     "glucocorticoids & antiepileptics",
                                     "estrogen"),
                      atc_codes = c('N06AB',
                                    'C10AA',
                                    'M01AH',
                                    'N05BA',
                                    'N02A',
                                    'N03A|N05AH02',
                                    'H02AB|N03AB02|N03AF01',
                                    'G03C'),
                      endpoints = c(vec2sep(bleed_endpoints),
                                   vec2sep(myopathy_endpoints),
                                   vec2sep(cvd_endpoints),
                                   vec2sep(bzd_endpoints),
                                   vec2sep(opioid_endpoints),
                                   vec2sep(agranulocytosis_endpoints),
                                   vec2sep(osteoporosis_endpoints),
                                   vec2sep(clot_endpoints)),
                              stringsAsFactors=FALSE)}

# read vnr_info, have to force VNRO as col_character() to prevent string"0012" from becoming integer "12"
vnr_info <- read_csv("./drugs/vnr_info_combined.csv", col_types = cols(VNRO = col_character()))

# save events endpoint_combination results to a list
events_list <- list()
tic()
for (i in 1:nrow(atc_endpoint_map)) {
    events_list[[atc_endpoint_map$endpoints[i]]] <- endpoint_combination(events, atc_endpoint_map[i,3])
}
toc()
events_save <- events


# combine exposure periods ignoring package size with default package size of 100
tolerance <- 14
ignore_psize <- TRUE
psize <- 100

exposures_list_ignore <- list()
for (i in 1:nrow(atc_endpoint_map)) {
    atc_codes <- data.frame(ATC_CODE = sep2vec(atc_endpoint_map[i,2]),
                        DAILY_DOSE = rep(1, length(sep2vec(atc_endpoint_map[i,2]))),
                       stringsAsFactors=FALSE)
    tic()
    print("Current drug_group: ")
    print(atc_endpoint_map[i,1])
    exposures_list_ignore[[atc_endpoint_map$drug_group[i]]] <- exposure_combination(purchases,
                                                                             atc_codes,
                                                                             vnr_info,
                                                                             tolerance,
                                                                             ignore_psize,
                                                                             psize)
    toc()
}

# combine exposure periods using package sizes
tolerance <- 14
ignore_psize <- FALSE

exposures_list <- list()
system.time(for (i in 1:nrow(atc_endpoint_map)) {
    atc_codes <- data.frame(ATC_CODE = sep2vec(atc_endpoint_map[i,2]),
                        DAILY_DOSE = rep(1, length(sep2vec(atc_endpoint_map[i,2]))),
                       stringsAsFactors=FALSE)
    tic()
    print("Current drug_group: ")
    print(atc_endpoint_map[i,1])
    exposures_list[[atc_endpoint_map$atc_codes[i]]] <- exposure_combination(purchases,
                                                                             atc_codes,
                                                                             vnr_info,
                                                                             tolerance,
                                                                             ignore_psize,
                                                                             psize)
    toc()
})