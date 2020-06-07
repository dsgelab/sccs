
library(tidyverse)
library(SCCS)


# endpoint_combination returns events included in vector "endpoints"
# events     : data frame of all recorded events (FINNGENID, EVENT_AGE, ENDPOINT)
# endpoints  : vector of FinnGen endpoint names to be studied
endpoint_combination <- function(events, endpoints) {
  
  # calculate observation start age
  # TODO: replace with correct years
  events$OBSERV_START_AGE <- events$FU_END_AGE - (2018-1995)
  
  # sort events by finngenid and event age
  events <- arrange(events, FINNGENID, EVENT_AGE)
  
  # filter events by "endpoints" vector
  events <- events[events$ENDPOINT %in% endpoints,]
  
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
  purchases <- purchases[purchases$ATC_CODE %in% atc_codes$ATC_CODE,]
  
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
  
  # create equally populated age groups, number specified in age_groups
  ageq <- floor(quantile(data$EVENT_AGE[duplicated(data$FINNGENID)==0],
                       seq(1/age_groups, 1-1/age_groups, 1/age_groups),names=F))
  
  # SCCS model fitting
  data$ref <- round(age_groups/2) # for passing reference age to function
  sccs.mod <- standardsccs(event~PURCHASE_AGE+relevel(age, ref=ref[1]),
                             indiv=FINNGENID, astart=OBSERV_START_AGE, aend=FU_END_AGE,
                             aevent=EVENT_AGE, adrug=PURCHASE_AGE, aedrug=EXPOSURE_END_AGE,
                             expogrp=exposure_periods, agegrp=ageq, dataformat="stack", data=data)
  
  return(sccs.mod)
}



# Manual work needed
# - diagnosis codes to vector "endpoints"
# - ATC codes and daily doses to dataframe "atc_codes"
# - choose exposure periods (standard is 0; if event-dependent exposure is suspected
# then add pre-exposure period e.g. c(-7, 0), not longer than tolerance)
# - choose number of age groups (could do hist(events$EVENT_AGE) first and play with
# the number of breaks so that most trends are accounted for but the division is not
# too fine to cause overfitting)
# - include_unexposed: usually TRUE, improves estimation of age effects
# - first_event: choose only first event or all events by individual
# - test differences in models with lrtsccs(model2, model1)


# read full event data, assign arguments and prepare relevant event data
events <- read.csv("mock_longitudinal_endpoint_data.txt", header=TRUE)
events$FU_END_AGE <- events$EVENT_AGE + 1 # this column was missing
endpoints <- c("E4_DIABETES", "COPD_COMORB")
events <- endpoint_combination(events, endpoints)

# read full exposure data, assign arguments and prepare relevant exposure data
purchases <- read.csv("mock_doac_purchases.csv", header=TRUE)
atc_codes <- data.frame(ATC_CODE = c("B01AF01", "B01AF02", "B01AE07"), DAILY_DOSE = c(1, 2, 2))
vnr_info <- read.csv("vnr_info_combined.csv", header=TRUE)
tolerance <- 14
ignore_psize <- TRUE
psize <- 100
exposures <- exposure_combination(purchases, atc_codes, vnr_info, tolerance, ignore_psize, psize)


# create single SCCS model
age_groups <- 20
hist(events$EVENT_AGE, breaks=quantile(events$EVENT_AGE, seq(0, 1, 1/age_groups)))
exposure_periods <- c(-7, 0)
include_unexposed <- TRUE
first_event <- FALSE
sccs.mod <- sccs_model(events, exposures, exposure_periods, age_groups,
           include_unexposed, first_event)
sccs.mod




# create array of SCCS models

# input list
inputs <- expand.grid(exposure_periods=list(0, c(-7, 0), c(-14, 0), c(0, 7)),
                      age_groups=c(10, 20, 30),
                      include_unexposed=c(TRUE,FALSE),
                      first_event=c(TRUE,FALSE))
inputs

# histograms with differing numbers of age groups
par(mfrow=c(1, length(unique(inputs$age_groups))))
for (n in unique(inputs$age_groups)) {
  hist(events$EVENT_AGE, breaks=quantile(events$EVENT_AGE, seq(0, 1, 1/n)), main=paste(n, "groups"))
}

# results_mod: list of models
# results_df: dataframe of inputs and relevant results
results_mod <- pmap(.l=inputs, .f=sccs_model, events=events, exposures=exposures)

results_df <- inputs
results_df[c("exp(coef)", "lower95", "upper95", "p", "se")] <- NA
for (row in 1:length(results_mod)) {
  i.exp <- length(results_df[[row,1]]) # index of relevant exposure group, assumed to be the last one
  results_df[row,5:9] <- with(results_mod[[row]], c(conf.int[i.exp,1], conf.int[i.exp,3],
                            conf.int[i.exp,4], coefficients[i.exp,5], coefficients[i.exp,3]))
}
results_df
