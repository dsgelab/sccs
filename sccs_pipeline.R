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
                              stringsAsFactors=FALSE)

# endpoint_combination returns events included in vector "endpoints"
# events     : data frame of all recorded events (FINNGENID, EVENT_AGE, ENDPOINT)
# endpoints  : vector of FinnGen endpoint names to be studied
endpoint_combination <- function(events, endpoints) {
    
    # modify | separated endpoints to vector
    endpoints <- sep2vec(endpoints)
    
    # calculate observation start age
    events$OBSERV_START_AGE <- events$FU_START_AGE

    # sort events by finngenid and event age
    events <- arrange(events, FINNGENID, EVENT_AGE)

    # filter events by "endpoints" vector
    events <- events[events$ENDPOINT %in% endpoints,]
    
    if (nrow(events)==0) {
        print("Zero rows for endpoint ", endpoint, ", aborting")
        return(NA)
    }

    return(events)
}



# exposure_combination returns exposures included in dataframe "atc_codes"
# purchases    : data frame of all purchases (FINNGENID, ATC_CODE, VNRO, PLKM, PURCHASE_AGE)
# atc_codes    : data frame of drugs to be studied (ATC_CODE, DAILY_DOSE)
# vnr_info     : data frame of package information (VNRO, PKOKO)
# tolerance    : time allowed between medication episodes that are to be considered one (in days) 
# ignore_psize : if TRUE, all purchases last for a default time, e.g. 100 days (use psize)
# psize        : default package size (in days)
exposure_combination <- function(purchases, atc_codes, vnr_info,
                                 tolerance=14, ignore_psize=FALSE, psize=100) {
    
    # filter purchases by relevant atc codes
    # returns all purchases with any of atc_codes$ATC_CODE at the beginning
    pattern <- paste("^", paste(atc_codes$ATC_CODE, collapse="|^"), sep="")
    purchases <- filter(purchases, grepl(pattern, ATC_CODE))
    
    if (nrow(purchases)==0) {
        print("Zero rows for ATCs ", pattern, ", aborting")
        return(NA)
    }

    # if given actual package sizes, use them
    if (!ignore_psize) {
    # combine purchases with package info
    purchases <- left_join(purchases, vnr_info, by="VNRO")

    # check if all VNRs were found
    if (any(is.na(purchases$PKOKO))) {
      print("VNRs not found:")
      print(unique(purchases[is.na(purchases$PKOKO),]$VNRO))
      stop()
    }
    print(head(purchases))
    # calculate exposure duration and exposure end age using doses per day
    purchases <- left_join(purchases, atc_codes, by="ATC_CODE")
    purchases$EXPOSURE_END_AGE <- purchases$PURCHASE_AGE +
      with(purchases, PLKM * PKOKO / 365 / DAILY_DOSE)

    # otherwise, assume default package size
    } else {
    purchases$EXPOSURE_END_AGE <- purchases$PURCHASE_AGE + psize/365
    }

    # sort purchases by finngenid and purchase age
    purchases <- arrange(purchases, FINNGENID, PURCHASE_AGE)

    # initialize episode number
    purchases$EPISODE <- NA
    purchases[1,]$EPISODE <- 1

    # loop over purchases and combine with previous purchase if within tolerance
    for (row in 2:nrow(purchases)) {
    # if same finngenid as previous:
    if (purchases[row,]$FINNGENID == purchases[row-1,]$FINNGENID) {
      # if the next purchase made within tolerance:
      if (purchases[row,]$PURCHASE_AGE < purchases[row-1,]$EXPOSURE_END_AGE + tolerance/365) {
        # keep episode number
        purchases[row,]$EPISODE <- purchases[row-1,]$EPISODE
        # if the next purchase made before previous has run out:
        if (purchases[row,]$PURCHASE_AGE < purchases[row-1,]$EXPOSURE_END_AGE) {
          # transfer excess exposure duration if actual package sizes are used
          if (!ignore_psize) {
            purchases[row,]$EXPOSURE_END_AGE <- purchases[row,]$EXPOSURE_END_AGE + 
              (purchases[row-1,]$EXPOSURE_END_AGE - purchases[row,]$PURCHASE_AGE)
          }
        }
        # next purchase (do not run the line below)
        next
      }
    }
    # if not the same finngenid or purchase not made within tolerance, increment episode number
    purchases[row,]$EPISODE <- purchases[row-1,]$EPISODE + 1
    }

    # leave first purchase of episode
    combined_purchases <- purchases %>% distinct(EPISODE, .keep_all = TRUE)

    # combined exposure end age <- end age of last exposure
    combined_purchases$EXPOSURE_END_AGE <- purchases %>%
    group_by(EPISODE) %>%
    summarize(end = max(EXPOSURE_END_AGE)) %>%
    pull(end)

    # plkm and pkoko no longer relevant
    combined_purchases <- subset(combined_purchases, select = c(FINNGENID, ATC_CODE, VNRO,
                                                    PURCHASE_AGE, EXPOSURE_END_AGE))

    return(combined_purchases)
}



# sccs_model returns the fitted SCCS model
# events            : data frame returned by endpoint_combination
# exposures         : data frame returned by exposure_definition
# exposure_periods  : vector of days to the start of exposure-related risk, e.g.
#                       0: no pre-exposure period, risk period starts at PURCHASE_AGE
#                       c(-7, 0): 7-day pre-exposure period
# age_groups        : number of age groups
# include_unexposed : TRUE: include individuals that have had event but no exposure
# first_event       : TRUE: select only first event per individual
sccs_model <- function(events, exposures, exposure_periods, age_groups,
                       include_unexposed, first_event) {
  
  # if first_event, filter subsequent events
  if (first_event) {
    events <- distinct(events, FINNGENID, .keep_all = TRUE)
  }
    print(head(exposures))
  
  # combine event and exposure data
  # if an individual can have multiple events, all exposure episodes are repeated for each event
  # if an individual has an event but no exposures, exposure-related columns will have NA values
  data <- left_join(events, exposures, by="FINNGENID")
  
  # if unexposed are not to be included, omit NA
  if (!include_unexposed) {
    data <- na.omit(data)
  }
  
  # convert years to integer days
  data <- data %>%
    mutate_at(vars(EVENT_AGE, FU_END_AGE, OBSERV_START_AGE, PURCHASE_AGE,
                   EXPOSURE_END_AGE), .funs = funs(round(. * 365)))
    
    #data$OBSERV_START_AGE <- ifelse(data$OBSERV_START_AGE<0, 0, data$OBSERV_START_AGE)
  
  # create equally populated age groups, number specified in age_groups
  ageq <- floor(quantile(data$EVENT_AGE[duplicated(data$FINNGENID)==0],
                       seq(1/age_groups, 1-1/age_groups, 1/age_groups),names=F))
  print(head(data))
    print(summary(data))
    print(names(data))
    print(object.size(data))
  # SCCS model fitting
  #ref1 <- round(age_groups/2) # for passing reference age to function
  sccs.mod <- standardsccs(event~PURCHASE_AGE+relevel(age, ref=5),
                             indiv=FINNGENID, astart=OBSERV_START_AGE, aend=FU_END_AGE,
                             aevent=EVENT_AGE, adrug=PURCHASE_AGE, aedrug=EXPOSURE_END_AGE,
                             expogrp=exposure_periods, agegrp=ageq, dataformat="stack", data=data)
  
  return(sccs.mod)
}