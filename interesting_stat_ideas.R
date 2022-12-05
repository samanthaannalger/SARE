# SARE Field Project Data Analysis
# P. Alexander Burnham
# 28 November 2022

# set directory:
setwd("~/Documents/GitHub/SARE")


# install libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(tidyr)
library(viridis)
library(car)
library(imputeTS)


# load in data
#ds <- read.csv("SARE_Field_database2022.csv", header = TRUE, stringsAsFactors = FALSE)
ds <- read.csv("SARE_field_database2022.csv", header = TRUE, stringsAsFactors = FALSE)
virus <- read.csv("DWV_SARE2021.csv", header = TRUE, stringsAsFactors = FALSE)

# colonies that were removed
d <- ds[grepl(ds$comments, pattern = "removed from", fixed = TRUE),]

# make sure ids are unique
unique_to_remove <- unique(d$lab_ID)

# pull rows out that match these values
ds = filter(ds, !(lab_ID %in% unique_to_remove))

# this is where we split by year
ds_2021<- ds[ds$year == 2021, ]
ds_2022 <- ds[ds$year == 2022, ]

# take only usefull columns
ds_2022_clean <- select(ds_2022, yard, treatment_grp, lab_ID, sampling_event, overwinter_success, frame_of_bees_FHA, FKB_percentile, UBO_assay_score, varroa_load_mites.100.bees, nosema_load_spores.bee)

# remove NA rows
ds_2022_clean <- ds_2022_clean[rowSums(is.na(ds_2022_clean)) != ncol(ds_2022_clean), ]

# fk and ubo binary variables created here
ds_2022_clean$FK_binary <- ifelse(ds_2022_clean$FKB_percentile >= 0.95, 1, 0) 
ds_2022_clean$UBO_binary <- ifelse(ds_2022_clean$UBO_assay_score >= 0.65, 1, 0) 


# summarize data set into one time point
summary_2022 <- ds_2022_clean %>% 
  group_by(yard, lab_ID, overwinter_success) %>% 
  summarise(
    
    varroa = max(varroa_load_mites.100.bees, na.rm=T), 
    nosema = max(nosema_load_spores.bee, na.rm=T), 
    FKA = sum(FK_binary, na.rm=T),
    UBO = sum(UBO_binary, na.rm=T),
    FOB = max(frame_of_bees_FHA, na.rm=T)
    
  ) 



# remove dead and missing colonies
summary_2022_clean <- summary_2022[!grepl("dead", summary_2022$overwinter_success),]

# remove column: overwinter success
summary_2022_clean$overwinter_success =  NULL

# rename sticker 8 to 8 and make column numeric
summary_2022_clean[summary_2022_clean$lab_ID == 56,]$FOB <- "8"
summary_2022_clean$FOB <- as.numeric(summary_2022_clean$FOB)

View(ds_2022)


