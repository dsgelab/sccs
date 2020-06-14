
purchases <- read.csv("/Users/jsjukara/Dropbox/FIMM/Sami Kulju/mock_doac_purchases.csv")

library(profvis)
library(tictoc)
atc_codes <- data.frame(ATC_CODE = c("B01AF01", "B01AF02", "B01AE07"), DAILY_DOSE = c(1, 2, 2))
vnr_info <- read.csv("/Users/jsjukara/Dropbox/FIMM/Sami Kulju/vnr_info_combined.csv", header=TRUE)
tolerance <- 14
ignore_psize <- TRUE
psize <- 100
library(tidyverse)

tic()
profvis({
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
  
  
  exposures <- exposure_combination(purchases, atc_codes, vnr_info, tolerance, ignore_psize, psize) })

toc()
exposures_save <- exposures

tic()
profvis({
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
    purchases$EPISODE[1] <- 1
    
    # loop over purchases and combine with previous purchase if within tolerance
    for (row in 2:nrow(purchases)) {
      # if same finngenid as previous:
      if (purchases$FINNGENID[row] == purchases$FINNGENID[row-1]) {
        # if the next purchase made within tolerance:
        if (purchases$PURCHASE_AGE[row] < purchases$EXPOSURE_END_AGE[row-1] + tolerance/365) {
          # keep episode number
          purchases$EPISODE[row] <- purchases$EPISODE[row-1]
          # if the next purchase made before previous has run out:
          if (purchases$PURCHASE_AGE[row] < purchases$EXPOSURE_END_AGE[row]) {
            # transfer excess exposure duration if actual package sizes are used
            if (!ignore_psize) {
              purchases$EXPOSURE_END_AGE[row] <- purchases$EXPOSURE_END_AGE[row] + 
                (purchases$EXPOSURE_END_AGE[row-1] - purchases$PURCHASE_AGE[row])
            }
          }
          # next purchase (do not run the line below)
          next
        }
      }
      # if not the same finngenid or purchase not made within tolerance, increment episode number
      purchases$EPISODE[row] <- purchases$EPISODE[row-1] + 1
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
  
  
  exposures <- exposure_combination(purchases, atc_codes, vnr_info, tolerance, ignore_psize, psize) })
toc()

identical(exposures_save,exposures)
all(exposures_save == exposures)

117.658/13.188
