
# Cornell Pesticide Data Analysis 

library(tidyverse)
library(dplyr)
library(ggplot2)

setwd("~/Documents/GitHub/SARE")

# Read in Tosi Datasets 
tosi_lethal <- read.csv("pesticide_data_to_merge/Tosi_lethal.csv", header = TRUE, stringsAsFactors = FALSE)
tosi_sublethal <- read.csv("pesticide_data_to_merge/Tosi_sublethal.csv", header = TRUE, stringsAsFactors = FALSE)

# Read in Description Dataset 
pest_Desc <- read.csv("pesticide_data_to_merge/pestDesc.csv", header = TRUE, stringsAsFactors = FALSE)


# Read in Cornell Dataset !NOTE! - find more general solution to white space/column headers
pest_Results <- read.csv("pesticide_data_to_merge/pesticide_results_2021.csv", header = TRUE, 
                         stringsAsFactors = FALSE, skip = 1)


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


########################################################################
# LIMIT FINDER FUNCTION
########################################################################
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
########################################################################


LS_df <- limit_finder(df = LS_df, search = "<LOQ", lookup = LS_lookup)
SS_df <- limit_finder(df = SS_df, search = "<LOQ", lookup = SS_lookup)
LS_df <- limit_finder(df = LS_df, search = ">ULOQ", lookup = LS_lookup)
SS_df <- limit_finder(df = SS_df, search = ">ULOQ", lookup = SS_lookup)

# append dataframes
pest_df <- rbind(LS_df, SS_df)

# create small scale/large scale column
pest_df$scale <- ifelse(pest_df$Mass..g. < 1, "small", "large")

str(pest_df)






################################################################################
# Cleaning LD50 Dataset -- Tosi Lethal
################################################################################

# looking to find the min LD50 value whether it be from contact or acute exposure types 
# put min value in its own column 

view(tosi_lethal)

# convert blank spaces to NA 
tosi_lethal[tosi_lethal == " "] <- NA
tosi_lethal[tosi_lethal == ""] <- NA

# adjusting future column names 
tosi_lethal["X.6"][tosi_lethal["X.6"] == "Min (ug/bee)"] <- "oral_acute_LD50_min"
tosi_lethal["X.13"][tosi_lethal["X.13"] == "LD50 1"] <- "oral_acute_LD50_1"
tosi_lethal["X.14"][tosi_lethal["X.14"] == "LD50 2"] <- "oral_acute_LD50_2"
tosi_lethal["X.15"][tosi_lethal["X.15"] == "LD50 3"] <- "oral_acute_LD50_3"
tosi_lethal["X.16"][tosi_lethal["X.16"] == "LD50 4"] <- "oral_acute_LD50_4"
tosi_lethal["X.17"][tosi_lethal["X.17"] == "LD50 5"] <- "oral_acute_LD50_5"

tosi_lethal["X.19"][tosi_lethal["X.19"] == "Min (ug/bee)"] <- "contact_acute_LD50_min"
tosi_lethal["X.26"][tosi_lethal["X.26"] == "LD50 1"] <- "contact_acute_LD50_1"
tosi_lethal["X.27"][tosi_lethal["X.27"] == "LD50 2"] <- "contact_acute_LD50_2"
tosi_lethal["X.28"][tosi_lethal["X.28"] == "LD50 3"] <- "contact_acute_LD50_3"

view(tosi_lethal)

# adjusting header of columns 
names(tosi_lethal) <- tosi_lethal[1,]
tosi_lethal <- tosi_lethal[-1,]

# renaming columns 
names(tosi_lethal)[names(tosi_lethal) == "Pesticide name"] <- "pesticide_name"


# find min LD50 value from both oral acute and contact acute values 
tosi_lethal$min_LD50_value <- min('oral_acute_LD50_min', 
                                   'oral_acute_LD50_1',
                                   'oral_acute_LD50_2',
                                   'oral_acute_LD50_3',
                                   'oral_acute_LD50_4',
                                   'oral_acute_LD50_5',
                                   'contact_acute_Ld50_min', 
                                   'contact_acute_LD50_1',
                                   'contact_acute_LD50_2',
                                   'contact_acute_LD50_3')

              # not printing min value 
                                   
view(tosi_lethal)


################################################################################
# Cleaning LOAEL Dataset -- Tosi Sublethal 
################################################################################

view(tosi_sublethal)

# renaming columns 
names(tosi_sublethal)[names(tosi_sublethal) == "ï..Pesticide.name"] <- "pesticide_name"
names(tosi_sublethal)[names(tosi_sublethal) == "LOAEL.Unit.measure"] <- "LOAEL_unit_measure"
names(tosi_sublethal)[names(tosi_sublethal) == "LOAEL.Lowest.sublethal.significant.dose.of.the.publication.per.effect.and.exposure.type..all.Unit.Measures."] <- "LOAEL"

# if LOAEL_unit_measure does not equal ppb, convert the values of LOAEL to ppb based on the unit measure of LOAEL_unit_measure.
## ex./ if LOAEL_unit_measure == ppm, then multiply LOAEL by 1000. (output in new column?)
## but if it is another unit, apply different conversion factor


# find the min LOAEL value among LOAELs for each pesticide (pesticide may occur more than once in rows)
## if possible between bee types (apis vs. anything else)

view(tosi_sublethal)
