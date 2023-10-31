# SARE Field Project Data Analysis
# P. Alexander Burnham
# 29 November 2021


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
virus2022 <- read.csv("DWV_SARE2022.csv", header = TRUE, stringsAsFactors = FALSE)

# colonies that were removed
d <- ds[grepl(ds$comments, pattern = "removed from", fixed = TRUE),]

# make sure ids are unique
unique_to_remove <- unique(d$lab_ID)

# pull rows out that match these values
ds = filter(ds, !(lab_ID %in% unique_to_remove))

# this is where we split by year
ds_2021<- ds[ds$year == 2021, ]
ds_2022 <- ds[ds$year == 2022, ]

# for 2022 create a column that denotes whether a colony was treated for varroa mites
# check sampling event number for future analyses
ds_2022$varroaTreated <- ifelse(ds_2022$sampling_event == 9 & ds_2022$varroa_load_mites.100.bees >= 2, "treated", "not_treated")

# removing colonies that were treated
ds_2022 <- ds_2022[ds_2022$varroaTreated == "not_treated",]

# merge in 2022 virus data and make log trans and prevalence data
ds_2022 <- merge(x = ds_2022,y = virus2022, all.x = TRUE)
ds_2022$NormGenomeCopy <- ifelse(is.na(ds_2022$NormGenomeCopy), 0, ds_2022$NormGenomeCopy)
ds_2022$virusPrev <- ifelse(ds_2022$NormGenomeCopy > 0, 1, 0) # making all negative IDs into 0s Chekc if there should be some NAs
ds_2022$logDWV <- log10(ds_2022$NormGenomeCopy + 1)

ds_2022$uboPrev <- ifelse(ds_2022$UBO_assay_score >= 0.6, 1, 0)


############################################################################
# 2022 Virus
############################################################################


# remove NAs from uba data set and make character variable
ds_2022UBO <- ds_2022[!is.na(ds_2022$uboPrev),]
ds_2022UBO$uboChar <- ifelse(ds_2022UBO$uboPrev==0, "UBO Neg.", "UBO pos.") 
no0 <- ds_2022UBO[ds_2022UBO$virusPrev == 1, ] # only positive virus bees

# ubo by virus load
ggplot(no0, aes(x=yard, y=NormGenomeCopy, color=uboChar)) + 
  geom_boxplot(size=1) +
  ylab("DWV (genome copies/bee)") + # y axis label
  xlab("Bee Yard") + # x axis label
  labs(color=NULL ) +
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("blue", "slategrey")) +
  scale_y_continuous(trans='log10')
  

vl <- lm(data = no0, logDWV ~ uboPrev * yard)
Anova(vl)



# summary of prevalence for virus by ubo and yard
prevSum <- ds_2022UBO %>% # operate on the dataframe (ds) and assign to new object (V)
  group_by(uboChar, yard) %>% # pick variables to group by
  summarise(
    mean = mean(virusPrev, na.rm=T), # mean
    N = length(virusPrev), # sample size
  ) 


# ub0 by virus prev
ggplot(prevSum, aes(x=yard, y=mean, fill=uboChar)) + 
  geom_bar(stat = "identity", position = "dodge") +
  ylab("DWV Prevalence") + # y axis label
  xlab("Bee Yard") + # x axis label
  labs(fill=NULL ) +
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_fill_manual(values = c("blue", "slategrey"))


vp <- glm(data = ds_2022UBO, virusPrev ~ uboPrev * yard, family = binomial(link="logit"))
Anova(vp)













############################################################################
# 2021 Analysis
############################################################################

#### VARROA ANALYSIS
# aggregate mite load by sampling event and yard
pltV2021 <- ds_2021 %>% # operate on the dataframe (ds) and assign to new object (V)
  group_by(sampling_event, yard) %>% # pick variables to group by
  summarise(
    
            mean = mean(varroa_load_mites.100.bees, na.rm=T), # mean
            SD = sd(varroa_load_mites.100.bees, na.rm=T), # standard dev.
            N = length(varroa_load_mites.100.bees), # sample size
            SE = SD/sqrt(N),                   # standard error
            MAX = max(varroa_load_mites.100.bees, na.rm=T)
            
            ) 


# Plot the time series data by group (yard in this case)
# the fist line of code calls in the data set and sets the variables
ggplot(data=pltV2021, aes(x=sampling_event, y=mean, group=yard, color=yard)) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Sampling Event") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  geom_line(size=1.5) + # create lines and set thickness 
  geom_point(size=4, shape=18) + # create points and set size and shape
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2) + # add standard errors
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Yard:") # color pallets option = A-H


## !!!How to group by month collected to show difference between varroa load progression of nucs vs full production over two years??
pltVtotal <- ds %>% # operate on the dataframe (ds) and assign to new object (Vtotal)
  pltVtotal$month <- format(as.Date(pltVtotal$date_collected, format="%m/%d/%Y"),"%m") %>%
  group_by(month, yard)  # pick variables to group by
summarise(
    
    mean = mean(varroa_load_mites.100.bees, na.rm=T), # mean
    SD = sd(varroa_load_mites.100.bees, na.rm=T), # standard dev.
    N = length(varroa_load_mites.100.bees), # sample size
    SE = SD/sqrt(N),                   # standard error
    MAX = max(varroa_load_mites.100.bees, na.rm=T)
    
  ) 



ggplot(data=pltVtotal, aes(x=sampling_event, y=mean, group=yard, color=yard)) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Sampling Event") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  geom_line(size=1.5) + # create lines and set thickness 
  geom_point(size=4, shape=18) + # create points and set size and shape
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2) + # add standard errors
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Yard:") # color pallets option = A-H




# statistical analysis 
# create the model
mod1 <- lmer(data=ds_2021, varroa_load_mites.100.bees ~ yard * sampling_event + (1|lab_ID))
summary(mod1) # look at the summary of the model
Anova(mod1) # check significance 



#### NOSEMA ANALYSIS
# aggregate nosema load by sampling event and yard
pltN <- ds_2021 %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(sampling_event, yard) %>% # pick variables to group by
  summarise(
    
    mean = mean(nosema_load_spores.bee, na.rm=T), # mean
    SD = sd(nosema_load_spores.bee, na.rm=T), # standard dev.
    N = length(nosema_load_spores.bee), # sample size
    SE = SD/sqrt(N),                   # standard error
    MAX = max(nosema_load_spores.bee, na.rm=T)
    
  ) 

# Plot the time series data by group (yard in this case)
# the fist line of code calls in the data set and sets the variables
ggplot(data=pltN, aes(x=sampling_event, y=mean, group=yard, color=yard)) +
  ylab("Nosema Load (spores/bee)") + # y axis label
  xlab("Sampling Event") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  geom_line(size=1.5) + # create lines and set thickness 
  geom_point(size=4, shape=18) + # create points and set size and shape
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2) + # add standard errors
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Yard:") # color pallets option = A-H

# statistical analysis 
# create the model
mod1 <- lmer(data=ds, nosema_load_spores.bee ~ yard * sampling_event + (1|lab_ID))
summary(mod1) # look at the summary of the model
Anova(mod1) # check significance 



#### FKA 2021 ANALYSIS
## here we are adding hygienic behavior to every instance of each

# subset of dataset with rows that have a LN2 test
FKB_2021 <- ds_2021[!is.na(ds$FKB_percentile),]


# select only columns we need
FKB_2021 <- select(FKB_2021, lab_ID, FKB_percentile)

# change column names
colnames(FKB_2021) <- c("lab_ID", "percent_hygienic")

# merge hygienic behavior back in
ds_2021_merged <- merge(x=ds_2021[!is.na(ds_2021$lab_ID),], y=FKB_2021[!is.na(FKB_2021$lab_ID),], by = "lab_ID", all.x = TRUE)



## Plot HB by varroa load
# Add regression lines
ggplot(ds_2021_merged, aes(x=percent_hygienic, y=varroa_load_mites.100.bees, 
               color=as.character(sampling_event))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))


# new plot that removes time points as a factor
ggplot(ds_2021_merged, aes(x=percent_hygienic, y=varroa_load_mites.100.bees)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) # size of the text and label ticks


# HB by varroa load model
mod2 <- lm(data=ds, varroa_load_mites.100.bees ~ FKB_percentile + sampling_event )
summary(mod2)
Anova(mod2)


# plot Hygienic behavior by yard

# not sure if needed here--->
# take only sampling event 2
# dst2 <- ds[ds$sampling_event==2,]

# Hygienic behavior by hive/yard, points number of hives, threshold dotted bar
ggplot(ds_2021_merged, aes(x=yard, y=percent_hygienic, color=yard)) + 
  geom_boxplot(size=1) +
  geom_text(aes(label=lab_ID), size=5) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  ylab("Percent Hygienic Behavior") + # y axis label
  xlab("Bee Yard") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H") +# color pallets option = A-H
  geom_hline(yintercept=.8, linetype="dashed",
             color = "red", size=1)


# set up the model
mod3 <- aov(data = ds_2021_merged, percent_hygienic ~ yard)
summary(mod3)


# print table in increasing order based on some variable
ds_2021_merged[order(ds_2021_merged$percent_hygienic, decreasing = TRUE),]



##### HONEY YIELD ANALYSIS
# take only sampling event 4
dst4 <- ds_2021_merged[ds_2021_merged$sampling_event==4,]

# plot Honey yield by hive/yard
ggplot(dst4, aes(x=yard, y=honey_removed, color=yard)) + 
  geom_boxplot(size=1) +
  geom_text(aes(label=lab_ID), size=5) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  ylab("Honey Yield (lbs.)") + # y axis label
  xlab("Bee Yard") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H") # color pallets option = A-H



#### FRAMES OF BEES ANALYSIS
# take only sampling event 1
dst1 <- ds_2021_merged[ds_2021_merged$sampling_event==1,]

# plot Honey yield by hive/yard
ggplot(dst1, aes(x=yard, y=frame_of_bees, color=yard)) + 
  geom_boxplot(size=1) +
  geom_text(aes(label=lab_ID), size=5) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  ylab("Frames of Bees") + # y axis label
  xlab("Bee Yard") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H") # color pallets option = A-H



# 2021 percent hygienic freeze kill
ds_2021_merged$FK_binary <- ifelse(ds_2021_merged$FKB_percentile >= 0.95, 1,0) #"hygienic", "nonhygienic")
mean(ds_2021_merged$FK_binary, na.rm=T) # get percentage of hygienic UBO




#########################################################################################################
# 2021 Survival Model

# create survival variable
ds_2021_merged$survival <- ifelse(ds_2021_merged$overwinter_success == "dead", 0, 
                                  ifelse(ds_2021_merged$overwinter_success == "alive", 1, NA))


survDF <- ds_2021_merged %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(field_ID) %>% # pick variables to group by
  summarise(
    
    varroa = mean(varroa_load_mites.100.bees, na.rm=T),
    nosema = mean(nosema_load_spores.bee, na.rm=T),
    hygienic = mean(FKB_percentile, na.rm=T),
    honey = mean(honey_removed, na.rm=T),
    weight = mean(october_weight, na.rm=T),
    surv = max(survival, na.rm=T),
    
)

survDF

max(ds_2021_merged$survival, na.rm=T)

table(ds_2021_merged$survival,ds_2021_merged$field_ID)







##############################################################################
# 2022 Data Analysis
##############################################################################

#### FKA 2022 ANALYSIS
## here we are adding hygienic behavior to every instance of each

# subset of dataset with rows that have a LN2 test
FKB_2022 <- ds_2022[!is.na(ds_2022$FKB_percentile),]
UBO_2022 <- ds_2022[!is.na(ds_2022$UBO_assay_score),]

# select only columns we need
FKB_2022 <- select(FKB_2022, lab_ID, FKB_percentile)
UBO_2022 <- select(UBO_2022, lab_ID, UBO_assay_score)

# change column names
colnames(FKB_2022) <- c("lab_ID", "percent_hygienic")
colnames(UBO_2022) <- c("lab_ID", "UBO_assay_score_merged")

# merge hygienic behavior back in
ds_2022 <- merge(x=ds_2022, y=FKB_2022, all.x = TRUE)
ds_2022 <- merge(x=ds_2022, y=UBO_2022, all.x = TRUE)

## Plot freeze kill by varroa load (continuous/ all varroa sampling events)
# Add regression lines
ggplot(ds_2022, aes(x=percent_hygienic, y=varroa_load_mites.100.bees, 
                    color=as.character(sampling_event))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("coral2", "grey", "deepskyblue2"), name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

## Plot freeze kill by varroa load (continuous/ only sampling event 9)
# Add regression lines
ggplot(August_Data, aes(x=percent_hygienic, y=varroa_load_mites.100.bees, 
                    color=as.character(sampling_event))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  ggtitle("Varroa Load by FKB Assay Percent") +
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("deepskyblue2"), name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

# new plot that removes time points as a factor
ggplot(ds_2022, aes(x=percent_hygienic, y=varroa_load_mites.100.bees)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  ggtitle("FKB Assay Percent by Varroa Load (composite)") +
  theme_minimal(base_size = 17) # size of the text and label ticks




# create binary variable for Freeze Kill
ds_2022$FK_binary <- ifelse(ds_2022$percent_hygienic >= 0.95, "hygienic", "non-hygienic")
mean(ds_2022$FK_binary, na.rm=T) # get percentage of hygienic FK

# create binary variable for UBO 
ds_2022$UBO_binary <- ifelse(ds_2022$UBO_assay_score >= 0.60, "hygienic", "non-hygienic") #"hygienic", "nonhygienic")
mean(ds_2022$UBO_binary, na.rm=T) # get percentage of hygienic UBO

# create binary variable for UBO for merged UBO dataset
ds_2022$UBO_binary_merged <- ifelse(ds_2022$UBO_assay_score_merged >= 0.60, "hygienic", "non-hygienic") #"hygienic", "nonhygienic")
mean(ds_2022$UBO_binary, na.rm=T) # get percentage of hygienic UBO



# differences in freeze kill between selection processes
ggplot(ds, aes(x=treatment_grp, y=FKB_percentile, color=treatment_grp, na.rm=T))+
  geom_boxplot(size=1) + 
  ylab("Percent Hygienic Behavior") + # y axis label
  xlab("Treatment Group") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("darkseagreen4","darkorange"), name="Treatment Group:")# color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  
mod <- aov(ds_2022$percent_hygienic ~ ds_2022$treatment_grp)
summary(mod)




# differences in UBO score between selction processes
ggplot(ds_2022, aes(x=treatment_grp, y=UBO_assay_score, color=treatment_grp)) +
  geom_boxplot(size=1) +
  ylab("Percent Hygienic Behavior") + # y axis label
  xlab("Treatment Group") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("darkseagreen4","darkorange"), name="Treatment Group:")# color pallets option = A-H



#geom_text(aes(label=lab_ID)) +
mod <- aov(ds_2022$UBO_assay_score ~ ds_2022$treatment_grp)
summary(mod)

# fisher exaxt tests (non parametric chi-square tests for selection method by 
# prortion hygienic)
fisher.test(y = ds_2022$FK_binary, x = ds_2022$treatment_grp)
fisher.test(y = ds_2022$UBO_binary, x = ds_2022$treatment_grp)


# differences in varroa load between selction processes over time
ds_2022 %>% # operate on the dataframe (ds) and assign to new object (V)
  group_by(sampling_event, treatment_grp) %>% # pick variables to group by
  summarise(
    
    mean = mean(varroa_load_mites.100.bees, na.rm=T), # mean
    SD = sd(varroa_load_mites.100.bees, na.rm=T), # standard dev.
    N = length(varroa_load_mites.100.bees), # sample size
    SE = SD/sqrt(N),                   # standard error
    MAX = max(varroa_load_mites.100.bees, na.rm=T)
    
  )  %>%

# Plot varroa loads by treatment group over time
ggplot(aes(x=sampling_event, y=mean, group=treatment_grp, color=treatment_grp, na.rm=T)) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Sampling Event") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  geom_line(size=1.5) + # create lines and set thickness 
  geom_point(size=4, shape=18) + # create points and set size and shape
  geom_errorbar(aes(ymin=mean-SE, ymax=mean+SE), width=.2) + # add standard errors
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Treatment Group:") # color pallets option = A-H


ds_2022_samp9 <- ds_2022[ds_2022$sampling_event == 9, ]
ds_2022_samp9$UBO_binary <- ifelse(ds_2022_samp9$UBO_assay_score >= 0.60, "hygienic", "non-hygienic") #"hygienic", "nonhygienic")
mean(ds_2022_samp9$UBO_binary, na.rm=T) # get percentage of hygienic UBO


##Bar chart using only late season varroa loads to compare treatment groups
ggplot(ds_2022_samp9, aes(x=treatment_grp, y=varroa_load_mites.100.bees, color=treatment_grp)) +
  geom_boxplot(size=1) +
  ylab("Varroa Load") + # y axis label
  xlab("Treatment Group") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("darkseagreen4","darkorange"), name="Treatment Group:")

# varroa load by group FHA vs NPQ (ANOVA)
mod <- aov(ds_2022$varroa_load_mites.100.bees ~ ds_2022$treatment_grp)
summary(mod)



# varroa load by UBO Binary


 # mutate(ds_2022, UBO_binary = ifelse(UBO_assay_score > .60, "Hygienic", "Non-hygienic"))

ggplot(data = ds_2022, mapping = aes(x=UBO_binary, y=varroa_load_mites.100.bees, color=UBO_binary, na.rm=T))+
 geom_boxplot(size=1) + 
 ylab("Varroa Load (mites/100 bees)") + # y axis label
 xlab("UBO Hygienic Behavior") + # x axis label
 theme_minimal(base_size = 17) + # size of the text and label ticks
 theme(legend.position = "top") + # place the legend at the top
 scale_color_manual(values = c("darkgrey","darkgreen"), name="Behavior:")# color pallets option = A-H

summary(aov(ds_2022$varroa_load_mites.100.bees~ds_2022$UBO_binary))



# select august course
August_Data = ds_2022[ds_2022$sampling_event == 9,] # UBO SE7

# varroa load by FK Binary August
ggplot(August_Data %>% filter(!is.na(FK_binary)), aes(x=FK_binary, y=varroa_load_mites.100.bees, color=FK_binary, ))+
  geom_boxplot(size=1) + 
  ylab("August Varroa Load (mites/ 100 bees)") + # y axis label
  xlab("FKB Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("darkseagreen4","darkorange"), name="Behavior:")# color pallets option = A-H

summary(aov(August_Data$varroa_load_mites.100.bees ~ August_Data$FK_binary))

# varroa load by UBO Binary August
ggplot(August_Data %>% filter(!is.na(FK_binary)), aes(x=UBO_binary_merged, y=varroa_load_mites.100.bees, color=UBO_binary_merged))+
  geom_boxplot(size=1) + 
  ylab("August Varroa Load (mites/ 100 bees)") + # y axis label
  xlab("UBO binary") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("darkseagreen4","darkorange"), name="Behavior:")# color pallets option = A-H
#geom_text(aes(label=lab_ID)) +

summary(aov(August_Data$varroa_load_mites.100.bees ~ August_Data$UBO_binary_merged))

hist(August_Data$varroa_load_mites.100.bees)


## Plot Varroa load by UBO score (continous)
# Add regression lines
ggplot(ds_2022, aes(x=UBO_assay_score, y=varroa_load_mites.100.bees, 
                    color=as.character(sampling_event))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))


# new plot that removes time points as a factor
ggplot(ds_2022, aes(x=UBO_assay_score, y=varroa_load_mites.100.bees)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) # size of the text and label ticks






# NEW ANALYSIS!!! nosema load by group FHA vs NPQ (ANOVA)
mod4 <- aov(ds_2022$nosema_load_spores.bee ~ ds_2022$treatment_grp)
summary(mod4)

# NEW ANALYSIS!!! virus load by group FHA vs NPQ (ANOVA)
##mod5 <- aov(ds_2022$virus_count ~ ds_2022$treatment_grp)
##summary(mod5)
## no 2022 data in virus count

# ALREADY COMPLETED!!! FKA (binary/continuous) by varroa load (omit treated hives)

# ALREADY COMPLETED!!! UBO (binary/continuous) by varroa load (omit treated hives)

# NEW ANALYSIS!!! FKA (binary/continuous) by nosema load (omit treated hives)
## Plot freeze kill by nosema load (continuous/ all nosema sampling events)
# Add regression lines


ggplot(ds_2022, aes(x=percent_hygienic, y=nosema_load_spores.bee, 
                        color=as.character(sampling_event))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Nosema Load (spores/bee)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

# nosema load by FK Binary (all nosema sampling events)
boxplot(ds_2022$nosema_load_spores.bee ~ ds_2022$FK_binary)
summary(aov(ds_2022$nosema_load_spores.bee ~ ds_2022$FK_binary))

# NEW ANALYSIS!!! UBO (binary/continuous) by nosema load (omit treated hives)
## Plot UBO by nosema load (continuous)
# Add regression lines
ggplot(ds_2022, aes(x=UBO_assay_score, y=nosema_load_spores.bee, 
                    color=as.character(sampling_event))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Nosema Load (spores/bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

# nosema load by UBO Binary
boxplot(ds_2022$nosema_load_spores.bee~ds_2022$UBO_binary)
summary(aov(ds_2022$nosema_load_spores.bee~ds_2022$UBO_binary))


# NEW ANALYSIS!!! FKA (binary/continuous) by virus load (omit treated hives)
## Plot FKA by virus load (continuous)
ggplot(ds_2022, aes(x=percent_hygienic, y=virus_count, 
                    color=as.character(sampling_event))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Virus Count") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

## virus load by FKA Binary (all nosema sampling events)

boxplot(ds_2022$virus_count ~ ds_2022$FK_binary)
summary(aov(ds_2022$virus_count ~ ds_2022$FK_binary))


# NEW ANALYSIS!!! UBO (binary/continuous) by virus load (omit treated hives)
## Plot UBO by virus load (continuous)
# Add regression lines
ggplot(ds_2022, aes(x=UBO_assay_score, y=virus_count, 
                    color=as.character(sampling_event))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Virus Count") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

## virus load by UBO (Binary)
boxplot(ds_2022$virus_count~ds_2022$UBO_binary)
summary(aov(ds_2022$virus_count~ds_2022$UBO_binary))


# NEW ANALYSIS!!! FKA by UBO categorical and continuous 
## Plot FKA by UBO (continuous)
# Add regression lines
ggplot(ds_2022, aes(x=UBO_assay_score, y=percent_hygienic, 
                    color=as.character(sampling_event))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=3) +
  ylab("FKB Assay Percent") + # y axis label
  xlab("UBO Assay Percent") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = c(5,5)) + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))



cor.test(ds_2022[ds_2022$sampling_event==7,]$UBO_assay_score_merged, ds_2022[ds_2022$sampling_event==7,]$percent_hygienic,  method="spearman", exact = F)



## FKA by UBO (Binary)
mosaicplot(ds_2022$FK_binary~ds_2022$UBO_binary)
summary(aov(ds_2022$FK_binary~ds_2022$UBO_binary))



## Plot UBO by varroa load (continuous)
# Add regression lines
ggplot(ds_2022, aes(x=UBO_assay_score, y=varroa_load_mites.100.bees, 
                    color=as.character(sampling_event))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))


# MODEL NOT WORKING-- unexpected symbol in mod2
# HB by varroa load model
mod2 <- lm(data=ds_2022, varroa_load_mites.100.bees ~ percent_hygienic + sampling_event)
summary(mod2)

Anova(mod2)


# plot Hygienic behavior by yard

# not sure if needed here--->
# take only sampling event 2
# dst2 <- ds[ds$sampling_event==2,]

# FKA by hive/yard, points number of hives, threshold dotted bar
ggplot(ds_2022, aes(x=yard, y=percent_hygienic, color=yard)) + 
  geom_boxplot(size=1) +
  geom_text(aes(label=lab_ID), size=5) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  ylab("Percent Hygienic Behavior") + # y axis label
  xlab("Bee Yard") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H") +# color pallets option = A-H
  geom_hline(yintercept=.95, linetype="dashed",
             color = "red", size=1)


# set up the model
mod3 <- aov(data = ds_2022, percent_hygienic ~ yard)
summary(mod3)


# print table in increasing order based on some variable
ds_2022[order(ds_2022$percent_hygienic, decreasing = TRUE),]


# UBO by hive/yard, points number of hives, threshold dotted bar
ggplot(ds_2022, aes(x=yard, y=UBO_assay_score, color=yard)) + 
  geom_boxplot(size=1) +
  geom_text(aes(label=lab_ID), size=5) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  ylab("Percent Hygienic Behavior") + # y axis label
  xlab("Bee Yard") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H") +# color pallets option = A-H
  geom_hline(yintercept=.6, linetype="dashed",
             color = "red", size=1)




















##################################################################################
### Code Stop (section models below for spring)
##################################################################################


##################################################################################
# selection model
# clean data:

# aggregate mite load by sampling event and yard
modelDF <- ds_2021 %>% # operate on the dataframe (ds) and assign to new object 
  group_by(lab_ID, yard) %>% # pick variables to group by
  summarise(
    varroa = mean(varroa_load_mites.100.bees, na.rm=T),
    nosema = mean(nosema_load_spores.bee, na.rm=T),
    hygienic = mean(FKB_percentile, na.rm=T),
    honey = mean(honey_removed, na.rm=T),
    weight = mean(october_weight, na.rm=T)
    
  ) 


modelDF <- merge(modelDF, virus, by = c("lab_ID"))

# impute missing values:
imputedDF <- na_mean(modelDF[,3:8])

# merge ds back together with ID vars
imputedDF <- cbind(imputedDF, modelDF[,1:2])

# rescale columns 0 to 1

# create function to rescale
range01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}



# call the function on each var
varroaScaled <- range01(imputedDF$varroa)
nosemaScaled <- range01(imputedDF$nosema)
hygienicScaled <- range01(imputedDF$hygienic)
honeyScaled <- range01(imputedDF$honey)
weightScaled <- range01(imputedDF$weight)
dwvScaled <- range01(imputedDF$NormGenomeCopy)

mn <- mean(sum(nosemaScaled), sum(weightScaled), sum(varroaScaled), sum(honeyScaled), sum(hygienicScaled), sum(dwvScaled))

varroaScaled <- (varroaScaled/sum(varroaScaled))*mn
nosemaScaled <- (nosemaScaled/sum(nosemaScaled))*mn
hygienicScaled <- (hygienicScaled/sum(hygienicScaled))*mn
honeyScaled <- (honeyScaled/sum(honeyScaled))*mn
weightScaled <- (weightScaled/sum(weightScaled))*mn
dwvScaled <- (dwvScaled/sum(dwvScaled))*mn

# create histogram data set of scaled vals
histogram <- data.frame(varroaScaled, nosemaScaled, hygienicScaled, honeyScaled, weightScaled, dwvScaled)
histogramLong <- gather(histogram, condition, measurement, varroaScaled, nosemaScaled, hygienicScaled, honeyScaled, weightScaled, dwvScaled)


# plot histogram of all scaled values
ggplot(histogramLong ,aes(x=measurement, fill=condition)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  ylab("Frequency") + # y axis label
  xlab("Parameter Value") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "right") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Parameter") # color pallets option = A-H



#Brood production= 3
#Varroa load= 2
#Nosema load= 1
#Virus load= 1
#Hygienic behavior= 5
#Honey production= 4
#October weight= 4
# 48, 51, 44, 63, 61


#c(varroa=12, nosema=5, hygienic=12, honey=10, weight=8, dwv=1)
# create fitness function
Fitness = (-3 * varroaScaled) +
  (-3 * nosemaScaled) +
  (3 * hygienicScaled) +
  (-3 * dwvScaled) +
  (3 * honeyScaled) +
  (3 * weightScaled)



# Fitness = (-16.8 * varroaScaled) + 
#   (18.92 * nosemaScaled) + 
#   (0.424 * hygienicScaled) + 
#   (-8.432 * dwvScaled) +
#   (-43.82 * honeyScaled) + 
#   (0.55 * weightScaled) 

# rescale 0 to 1
imputedDF$Fitness <- range01(Fitness)
death <- select(ds[ds$sampling_event==6,], overwinter_success, lab_ID)

imputedDF <- merge(imputedDF, death)
tp6 <- read.csv("TP6Data.csv")

imputedDF <- merge(imputedDF, tp6)

# print the sorted data set
orderedDF <- imputedDF[order(imputedDF$Fitness, decreasing = TRUE),]  
orderedDF



#write.csv(orderedDF[orderedDF$overwinter_success=="alive",], "modelselectionSARE22.csv")


pcaDS <- orderedDF[!orderedDF$lab_ID==34,]
pcaDS <- select(pcaDS, overwinter_success, varroa, nosema, hygienic, honey, weight, NormGenomeCopy)


pcaDS$overwinter_success


mod <- aov(pcaDS$honey~pcaDS$overwinter_success)
summary(mod)

boxplot(pcaDS$honey~pcaDS$overwinter_success)


set.seed(123)
ind <- sample(2, nrow(pcaDS),
              replace = TRUE,
              prob = c(0.6, 0.4))

training <- pcaDS[ind==1,]
testing <- pcaDS[ind==2,]




linear <- lda(overwinter_success~., training)
linear$scaling*1000000000000

p1 <- predict(linear, training)$class
tab <- table(Predicted = p1, Actual = training$overwinter_success)
tab

sum(diag(tab))/sum(tab)

p2 <- predict(linear, testing)$class
tab1 <- table(Predicted = p2, Actual = testing$overwinter_success)
tab1
sum(diag(tab1))/sum(tab1)






plot(sort(imputedDF$Fitness))

hist(imputedDF$Fitness, breaks = 10)



# constructing a =n LDA

library(klaR)
library(psych)
library(MASS)
library(ggord)
library(devtools)


# ##########################################################################################
# # amounts in each reaction
# amount <- rep(c(10^-2,10^-3, 10^-4, 10^-5, 10^-6, 10^-7, 10^-8, 10^-9), each = 5) * 4.11*10^10
# 
# # read in standard curve data
# ct3 <- read.csv("Alger-SA-17238.csv")
# x <- split(ct3, ct3$Target.Name)
# ACTIN <- x$ACTIN$CT
# DWV <- x$DWV$CT
# DF <- data.frame(DWV, ACTIN)
# 
# # get standards
# res <- pcr_assess(DF,
#                     amount = amount,
#                     method = 'standard_curve',
#                   plot=T)
# 
# 
# # get efficiency
# pcr_assess(DF,
#                  amount = amount,
#                  reference_gene = 'ACTIN',
#                  method = 'efficiency')
# ##########################################################################################
# 
# 
#source("BurnhamFunctionsSARE.R")
# 
#virusData <- read.csv("virusData2022.csv")
#dilution <- read.csv("RNAdilutions2022.csv")
# 

# Program Body:
# preliminary cleaning -> removes duplicate rows and control data
#TempVarClean <- PrelimClean(data = virusData)

# merge data sets to inlude dilution data for Normalization:
#TempVarClean <- merge(TempVarClean, dilution, by = "ID", all.x = TRUE)

# normalize data set Viral Load
#TempVarClean <- VirusNorm(number_bees = 50, data = TempVarClean)

# normalize viral laod by actin
#TempVarClean <- actinNormal(data = TempVarClean)

# make binary variable and use threashld of ct for limit of detection: 
#TempVarClean <- CT_Threash(data = TempVarClean)

#finalVirusDF <- TempVarClean[TempVarClean$target_name=="DWV",]
#DWV_SARE2022 <- select(finalVirusDF, ID, NormGenomeCopy) 

#write.csv(DWV_SARE2022, "DWV_SARE2022.csv")



































##################################################################################
###-------------------------------2022 Spring pics-----------------------------###
##################################################################################
# selection model


# clean data:
ds_2022_clean <- select(ds_2022, varroa_load_mites.100.bees, nosema_load_spores.bee, FKB_percentile, UBO_assay_score, NormGenomeCopy, yard, lab_ID)
ds_2022_clean <- ds_2022_clean[!is.na(ds_2022_clean$lab_ID), ]

# aggregate mite load by sampling event and yard
modelDF <- ds_2022_clean %>% # operate on the dataframe (ds) and assign to new object 
  group_by(lab_ID, yard) %>% # pick variables to group by
  summarise(
    varroa = mean(varroa_load_mites.100.bees, na.rm=T),
    nosema = mean(nosema_load_spores.bee, na.rm=T),
    fkb = mean(FKB_percentile, na.rm=T),
    ubo = mean(UBO_assay_score, na.rm=T),
    dwv = mean(NormGenomeCopy, na.rm=T),
  ) 



# impute missing values:
imputedDF <- na_mean(modelDF[,3:7])

# merge ds back together with ID vars
imputedDF <- cbind(imputedDF, modelDF[,1:2])

# create function to rescale
range01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}


# call the function on each var
varroaScaled <- range01(imputedDF$varroa)
nosemaScaled <- range01(imputedDF$nosema)
fkbScaled <- range01(imputedDF$fkb)
uboScaled <- range01(imputedDF$ubo)
dwvScaled <- range01(imputedDF$dwv)


mn <- mean(sum(nosemaScaled), sum(fkbScaled), sum(varroaScaled), sum(uboScaled), sum(dwvScaled))

varroaScaled <- (varroaScaled/sum(varroaScaled))*mn
nosemaScaled <- (nosemaScaled/sum(nosemaScaled))*mn
fkbScaled <- (fkbScaled/sum(fkbScaled))*mn
uboScaled <- (uboScaled/sum(uboScaled))*mn
dwvScaled <- (dwvScaled/sum(dwvScaled))*mn

# create histogram data set of scaled vals
histogram <- data.frame(varroaScaled, nosemaScaled, fkbScaled, uboScaled, dwvScaled)
histogramLong <- gather(histogram, condition, measurement, varroaScaled, nosemaScaled, fkbScaled, uboScaled, dwvScaled)


# plot histogram of all scaled values
ggplot(histogramLong ,aes(x=measurement, fill=condition)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  ylab("Frequency") + # y axis label
  xlab("Parameter Value") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "right") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Parameter") # color pallets option = A-H



#Brood production= 3
#Varroa load= 2
#Nosema load= 1
#Virus load= 1
#Hygienic behavior= 5
#Honey production= 4
#October weight= 4
# 48, 51, 44, 63, 61


#c(varroa=12, nosema=5, hygienic=12, honey=10, weight=8, dwv=1)
# create fitness function
Fitness = (-3 * varroaScaled) +
  (-3 * nosemaScaled) +
  (6 * uboScaled) +
  (-3 * dwvScaled) +
  (3 * fkbScaled)



# Fitness = (-16.8 * varroaScaled) + 
#   (18.92 * nosemaScaled) + 
#   (0.424 * hygienicScaled) + 
#   (-8.432 * dwvScaled) +
#   (-43.82 * honeyScaled) + 
#   (0.55 * weightScaled) 

# rescale 0 to 1
imputedDF$Fitness <- range01(Fitness)



# print the sorted data set
orderedDF <- imputedDF[order(imputedDF$Fitness, decreasing = TRUE),]  
orderedDF



write.csv(orderedDF, "modelselectionSARE23.csv")


pcaDS <- orderedDF[!orderedDF$lab_ID==34,]
pcaDS <- select(pcaDS, overwinter_success, varroa, nosema, hygienic, honey, weight, NormGenomeCopy)




mod <- aov(pcaDS$honey~pcaDS$overwinter_success)
summary(mod)

boxplot(pcaDS$honey~pcaDS$overwinter_success)


set.seed(123)
ind <- sample(2, nrow(pcaDS),
              replace = TRUE,
              prob = c(0.6, 0.4))

training <- pcaDS[ind==1,]
testing <- pcaDS[ind==2,]




linear <- lda(overwinter_success~., training)
linear$scaling*1000000000000

p1 <- predict(linear, training)$class
tab <- table(Predicted = p1, Actual = training$overwinter_success)
tab

sum(diag(tab))/sum(tab)

p2 <- predict(linear, testing)$class
tab1 <- table(Predicted = p2, Actual = testing$overwinter_success)
tab1
sum(diag(tab1))/sum(tab1)






plot(sort(imputedDF$Fitness))

hist(imputedDF$Fitness, breaks = 10)



