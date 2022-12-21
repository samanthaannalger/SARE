
# Cornell Pesticide Data Analysis 

setwd("~/GitHub/Pesticide")

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

# Merge Cornell Dataset and Description Dataset 
#cornell_pesticide_data <- merge(cornell_pesticide_data, pest_Desc, by = pesticide_name)



# Read in Look-Up Table
look_up_ds <- read.csv("pestLookUp.csv", header = TRUE, stringsAsFactors = FALSE)

# Read in 2021_Data_Formatted
formatted_ds <- read.csv("pestFormatted.csv", header = TRUE, stringsAsFactors = FALSE)
