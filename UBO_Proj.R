# UBO Bee Breeder Project Data Analysis
# S. Miller & C. McKay
# 30 November 2022


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
ds <- read.csv("UBO_Data_2022.csv", header = TRUE, stringsAsFactors = FALSE)

#################################################################################
# UBO Analysis
#################################################################################

# UBO cont and binary 
# create binary variable for UBO 
ds$UBO_binary <- ifelse(ds$assay_score >= 0.6, 1, 0) #"hygienic", "nonhygienic")
mean(ds$UBO_binary, na.rm=T) # get percentage of hygienic UBO


## UBO by June Varroa loads continuous 
# Add regression lines
ggplot(ds, aes(x=assay_score, y=june_varroa_load_mites.100.bees, 
                    color=as.character(location))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("June Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

## UBO by June varroa load binary 
boxplot(ds$june_varroa_load_mites.100.bees~ds$UBO_binary)
summary(aov(ds$june_varroa_load_mites.100.bees~ds$UBO_binary))


## UBO by August Varroa loads continuous
# Add regression lines
ggplot(ds, aes(x=assay_score, y=august_varroa_load_mites.100.bees, 
               color=as.character(location))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("August Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))


## UBO by August Varroa loads binary
boxplot(ds$august_varroa_load_mites.100.bees~ds$UBO_binary)
summary(aov(ds$august_varroa_load_mites.100.bees~ds$UBO_binary))


## UBO June composite continuous 
# Add regression lines
ggplot(ds, aes(x=assay_score, y=june_varroa_load_mites.100.bees)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("June Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

## UBO August composite continuous 
# Add regression lines
ggplot(ds, aes(x=assay_score, y=august_varroa_load_mites.100.bees)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("June Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

#################################################################################
# Nosema Analysis
#################################################################################

## UBO by June Nosema load continuous 
# Add regression lines
ggplot(ds, aes(x=assay_score, y=june_nosema_load_spores.bee, 
               color=as.character(location))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("June Nosema Load (spores/bee)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

## UBO by June Nosema load binary 
boxplot(ds$june_nosema_load_spores.bee~ds$UBO_binary)
summary(aov(ds$june_nosema_load_spores.bee~ds$UBO_binary))


## UBO by August Nosema load continuous 
# Add regression lines
ggplot(ds, aes(x=assay_score, y=august_nosema_load_spores.bee, 
               color=as.character(location))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("August Nosema Load (spores/bee)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

## UBO by June Nosema load binary 
boxplot(ds$august_nosema_load_spores.bee~ds$UBO_binary)
summary(aov(ds$august_nosema_load_spores.bee~ds$UBO_binary))



# UBO by frames of brood 
ggplot(ds, aes(x=frames_of_brood, y=assay_score, 
               color=as.character(location))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("UBO Assay Score)") + # y axis label
  xlab("Colony Size") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))
