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


# redefine atc_endpoint_map with negative controls
atc_endpoint_map_with_neg <- data.frame(drug_group = c('ssri_bleed',
                                     'statin_myopathy',
                                     "cox2i_cvd",
                                     "benzo_falls",
                                     "opioid_constipation",
                                     "antiepileptics_clozapine_agranulocytosis",
                                     "glucocorticoids_antiepileptics_osteoporosis",
                                     "estrogen_dvtpe",
                                     "statin_bleed",
                                     "statin_constipation",
                                     "statin_dvtpe",
                                     "benzo_myopathy",
                                     "benzo_cvd",
                                     "benzo_constipation",
                                     "benzo_dvtpe",
                                     "benzo_osteoporosis",
                                     "estrogen_constipation",
                                     "estrogen_falls"),
                      atc_codes = c('N06AB',
                                    'C10AA',
                                    'M01AH',
                                    'N05BA',
                                    'N02A',
                                    'N03A|N05AH02',
                                    'H02AB|N03AB02|N03AF01',
                                    'G03C',
                                    'N06AB',
                                    'N06AB',
                                    'N06AB',
                                    'N05BA',
                                    'N05BA',
                                    'N05BA',
                                    'N05BA',
                                    'N05BA',
                                    'G03C',
                                    'G03C'),
                      endpoints = c(vec2sep(bleed_endpoints),
                                   vec2sep(myopathy_endpoints),
                                   vec2sep(cvd_endpoints),
                                   vec2sep(bzd_endpoints),
                                   vec2sep(opioid_endpoints),
                                   vec2sep(agranulocytosis_endpoints),
                                   vec2sep(osteoporosis_endpoints),
                                   vec2sep(clot_endpoints),
                                   vec2sep(bleed_endpoints),
                                   vec2sep(opioid_endpoints),
                                   vec2sep(clot_endpoints),
                                   vec2sep(myopathy_endpoints),
                                   vec2sep(cvd_endpoints),
                                   vec2sep(opioid_endpoints),
                                   vec2sep(clot_endpoints),
                                   vec2sep(osteoporosis_endpoints),
                                   vec2sep(opioid_endpoints),
                                   vec2sep(bzd_endpoints)),
                              stringsAsFactors=FALSE)


# input list
inputs <- expand.grid(exposure_periods=list(0, c(-14, 0)),
                      age_groups=20,
                      include_unexposed=c(FALSE),
                      first_event=c(TRUE,FALSE))

# initiate parallelization
library(furrr)
plan(multiprocess)

analyses_list_ignore <- list()
results_list_ignore <- list()

for (i in 1:nrow(atc_endpoint_map_with_neg)) {
    events <- events_list[[atc_endpoint_map_with_neg$endpoints[i]]]
    exposures <- exposures_list_ignore[[atc_endpoint_map_with_neg$atc_codes[i]]] 
    
    coefficients_list <- list()
    
    
    # results_mod: list of models
    # results_df: dataframe of inputs and relevant results
    # use future_pmap() to parallelize the analyses for different inputs
    results_mod <- future_pmap(.l=inputs, .f=sccs_model, events=events, exposures=exposures)

    results_df <- inputs
    results_df[c("exp(coef)", "lower95", "upper95", "p", "se")] <- NA
    for (row in 1:length(results_mod)) {
        i.exp <- length(results_df[[row,1]]) # index of relevant exposure group, assumed to be the last one
        results_df[row,5:9] <- with(results_mod[[row]], c(conf.int[i.exp,1], conf.int[i.exp,3],
                                conf.int[i.exp,4], coefficients[i.exp,5], coefficients[i.exp,3]))
        
        # get all coefficients for one analysis with multiple inputs
        coefficients_list[[paste(inputs[i,], collapse="_")]] <- results_mod[[row]]$coefficients
        
    }
    # save list of coefficients to a list of analyses
    analyses_list_ignore[[atc_endpoint_map_with_neg$drug_group[i]]] <- coefficients_list
    
    # save primary coefficients
    results_list_ignore[[atc_endpoint_map_with_neg$drug_group[i]]] <- results_df
}
toc()

# initialize list that will contain one list per analysis, which will contain data frames of the coefficients
results_list <- list()
analyses_list <- list()

for (i in 1:nrow(atc_endpoint_map_with_neg)) {
    events <- events_list[[atc_endpoint_map_with_neg$endpoints[i]]]
    exposures <- exposures_list[[atc_endpoint_map_with_neg$atc_codes[i]]] 
    
    coefficients_list <- list()
    
    
    # results_mod: list of models
    # results_df: dataframe of inputs and relevant results
    # use future_pmap() to parallelize the analyses for different inputs
    results_mod <- future_pmap(.l=inputs, .f=sccs_model, events=events, exposures=exposures)

    results_df <- inputs
    results_df[c("exp(coef)", "lower95", "upper95", "p", "se")] <- NA
    for (row in 1:length(results_mod)) {
        i.exp <- length(results_df[[row,1]]) # index of relevant exposure group, assumed to be the last one
        results_df[row,5:9] <- with(results_mod[[row]], c(conf.int[i.exp,1], conf.int[i.exp,3],
                                conf.int[i.exp,4], coefficients[i.exp,5], coefficients[i.exp,3]))
        
        # get all coefficients for one analysis with multiple inputs
        coefficients_list[[paste(inputs[i,], collapse="_")]] <- results_mod[[row]]$coefficients
        
    }
    # save list of coefficients to a list of analyses
    analyses_list[[atc_endpoint_map_with_neg$drug_group[i]]] <- coefficients_list
    
    # save 
    results_list[[atc_endpoint_map_with_neg$drug_group[i]]] <- results_df
    print(paste("Completed: ", round(i/nrow(atc_endpoint_map_with_neg),2)*100, "%"), sep="")
}

toc()

