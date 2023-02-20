
# Cornell Pesticide Data Analysis 

library(tidyverse)
library(dplyr)
library(ggplot2)

setwd("~/Documents/GitHub/SARE")


## Reading in datasets 

# Read in Tosi Datasets 
tosi_lethal <- read.csv("pesticide_data_to_merge/Tosi_lethal.csv", header = TRUE, stringsAsFactors = FALSE, skip = 1)

tosi_sublethal <- read.csv("pesticide_data_to_merge/Tosi_sublethal.csv", header = TRUE, stringsAsFactors = FALSE)

# Read in Description Dataset 
pest_Desc <- read.csv("pesticide_data_to_merge/pestDesc.csv", header = TRUE, stringsAsFactors = FALSE)


# Read in Cornell Results Dataset !NOTE! - find more general solution to white space/column headers
pest_Results <- read.csv("pesticide_data_to_merge/pesticide_results_2021.csv", header = TRUE, 
                         stringsAsFactors = FALSE, skip = 1)


###############################################################################
# Cleaning Cornell Results Dataset
###############################################################################

# replace n.d. with NA
pest_Results[pest_Results == "n.d."] <- NA

# find our cut points
LS_row <- which(pest_Results == "Large Scale", arr.ind=TRUE)[1]
SS_row <- which(pest_Results == "Small Scale", arr.ind=TRUE)[1]
SS_filename_row <- which(pest_Results == "File Name", arr.ind=TRUE)[1]
bottom_row <- which(pest_Results == "Results are in ppb.", arr.ind=TRUE)[1] #find regex solution Results are in *

# cut out our data frames
LS_df <- pest_Results[1:(LS_row-1),]
SS_df <- pest_Results[(SS_filename_row+1):(SS_row-1),]

# cut out look up tables
LS_lookup <- pest_Results[LS_row:(LS_row+3),]
SS_lookup <- pest_Results[SS_row:(SS_row+3),]


#######################################################
# LIMIT FINDER FUNCTION
#######################################################
limit_finder <- function(df, search, lookup, scale){
  
  # find where samples say <loq
  loqVals <- data.frame(which(df == search, arr.ind=TRUE))
  
  if(length(loqVals$row>0)){
    
    # pull out mass and lod and do out the division
    mass <- as.numeric(df$Mass..g.[loqVals$row])
    scaleNum <- ifelse(search==">ULOQ", 3, 1) # convert scale into row index
    print(scaleNum)
    lod <- as.numeric(lookup[scaleNum,loqVals$col])
    
    results <- lod/mass
    
    # assign results to index where loq was found
    for(i in 1:length(results)){
      df[loqVals$row[i], loqVals$col[i]] <- results[i]
    }
  }
  return(df)
}
#######################################################

LS_df <- limit_finder(df = LS_df, search = "<LOQ", lookup = LS_lookup)
SS_df <- limit_finder(df = SS_df, search = "<LOQ", lookup = SS_lookup)
LS_df <- limit_finder(df = LS_df, search = ">ULOQ", lookup = LS_lookup)
SS_df <- limit_finder(df = SS_df, search = ">ULOQ", lookup = SS_lookup)

# append dataframes
pest_df <- rbind(LS_df, SS_df)

# create small scale/large scale column
pest_df$scale <- ifelse(pest_df$Mass..g. < 1, "small", "large")

str(pest_df)

view(pest_df)


################################################################################
# Cleaning LD50 Dataset -- Tosi Lethal
################################################################################
tosi_lethal
# looking to find the min LD50 value whether it be from contact or acute exposure types 
# put min value in its own column 

view(tosi_lethal)

# convert blank spaces to NA 
tosi_lethal[tosi_lethal == " "] <- NA
tosi_lethal[tosi_lethal == ""] <- NA


# changing column names 
colnames(tosi_lethal) # original column names

tosi_lethal_colnames <- c("pesticide_name", "other_names","cas", "pesticide_type", "MoA_short", "MoA_classification_site_target", "oral_LD50_geometricmean_ugbee", "oral_source_num","oral_LD50_min", "oralQ1", "oralQ2_median", "oralQ3", "oral_LD50_max", "oral_range", "oral_source_name", "oral_LD50_1", "oral_LD50_2", "oral_LD50_3", "oral_LD50_4", "oral_LD50_5", "contact_LD50_geometricmean_ugbee","contact_source_num","contact_LD50_min", "contactQ1", "contactQ2_median", "contactQ3", "contact_LD50_max","contact_range", "contact_source_name", "contact_LD50_1", "contact_LD50_2","contact_LD50_3")
  
colnames(tosi_lethal) <- tosi_lethal_colnames

#colnames(tosi_lethal) # verify new column names

# finding minLD50 value
tosi_lethal %>% rowwise() %>% mutate(min_LD50_value = min(oral_LD50_min, oral_LD50_1, oral_LD50_2, oral_LD50_3, oral_LD50_4, oral_LD50_5, contact_LD50_min, contact_LD50_1, contact_LD50_2, contact_LD50_3))
# !NOTE! not sure if this is working 

# view(tosi_lethal)


################################################################################
# Cleaning LOAEL Dataset -- Tosi Sublethal 
################################################################################

view(tosi_sublethal)

colnames(tosi_sublethal)

tosi_sublethal_colnames <- c("pesticide_name", "cas", "pesticide_type", "MoA_short", "MoA_classification_site_taret", "survey_inclusion_name", "screened_in_survey", "num_survery_screenings", "oral_LD50_geometricmean_ugbee", "oral_source_name", "contact_LD50_geometricmean_ugbee", "contact_source_name", "LOAEL_allunits", "LOAEL_unit_measure", "LOAEL_ug/bee/day", "LOAEL_category_ug/bee", "SubTR_LOAEL/LD50", "SubTR_category", "exposure_type_oral_v_contact", "exposure_type_acute_v_chronic", "exposure_duration_h", "time_after_exposure_of_significant_effect_h", "feedtype_main_category", "feedtype_subcategory", "feedtype_concentration", "bee_type", "bee_genus", "bee_species", "bee_species_details", "sublethal_effect_main_category", "sublethal_effect_subcategory", "sublethal_effect_details", "original_ref", "ref_year", "review_ref")

colnames(tosi_sublethal) <- tosi_sublethal_colnames

# colnames(tosi_sublethal) # verify new column names

tosi_sublethal[tosi_sublethal == " "] <- NA
tosi_sublethal[tosi_sublethal == ""] <- NA

# if LOAEL_unit_measure does not equal ppb, convert the values of LOAEL to ppb based on the unit measure of LOAEL_unit_measure.
## ex./ if LOAEL_unit_measure == ppm, then multiply LOAEL by 1000. (output in new column?)
## but if it is another unit, apply different conversion factor

tosi_sublethal$LOAEL_unit_measure <- as.character(tosi_sublethal$LOAEL_unit_measure)

str(tosi_sublethal$LOAEL_unit_measure)
which(table(tosi_sublethal$LOAEL_unit_measure)>=1)

# could unit measures be put into a function for conversion to ppb? 
tosi_sublethal_unit_measures <- c("µg/bee", "µM", "g/bee/week", "g/ha", "g/hive", "g/hm-2", "gals/acre", "μg", "μg/bee", 
                                  "μg/bee/day", "μg/larva", "μL", "μL/bee", "μM", "kg/ha", "MFR", "mL/bee", "mL/colony", 
                                  "mM", "mm3 /bee", "ng/L", "ng/ml", "nM", "nmol/bee", "nmol/day/bee", "ppb", "ppm", "unclear")

# tosi_sublethal %>% 
 #  select(LOAEL_unit_measure, LOAEL_allunits) %>% 
#  mutate(LOAEL_calculated = LOAEL_allunits * 1000)
  
#  if (LOAEL_unit_measure == "ppb") {
#    mutate(LOAEL_calculated == LOAEL_allunits * 1000)
#    } else {(LOAEL_calculated == "")
#    }


tosi_sublethal$LOAEL_calculated <- 
  if(LOAEL_unit_measure = "ppm") {
    mutate(tosi_sublethal$LOAEL_allunits * 1000)
  } 
elseif (LOAEL_unit_measure = "other unit measure") {
  mutate(tosi_sublethal$LOAEL_allunits * x conversion factor)
}
elseif (LOAEL_unit_measure = "ppb"){#no change
}
view(tosi_sublethal)

# find the min LOAEL value among LOAELs for each pesticide (pesticide may occur more than once in rows)
## if possible between bee types (apis vs. anything else)

