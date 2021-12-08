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
ds <- read.csv("SARE_Field_database.csv", header = TRUE, stringsAsFactors = FALSE)


# colonies that were removed
d <- ds[grepl(ds$comments, pattern = "removed from", fixed = TRUE),]


# make sure ids are unique
unique_to_remove <- unique(d$lab_ID)

# pull rows out that match these values
ds = filter(ds, !(lab_ID %in% unique_to_remove))


# aggregate mite load by sampling event and yard
pltDF <- ds %>% # operate on the dataframe (ds) and assign to new object (pltDF)
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
ggplot(data=pltDF, aes(x=sampling_event, y=mean, group=yard, color=yard)) +
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
mod1 <- lmer(data=ds, varroa_load_mites.100.bees ~ yard * sampling_event + (1|lab_ID))

summary(mod1) # look at the summary of the model

Anova(mod1) # check significance 



#### here we are adding hygienic behavior to every instnace of each I

# subset of dataset with rows that have a LN2 test
testedDS <- ds[!is.na(ds$HB_percentile),]
# select only columns we need
testedDS <- select(testedDS, lab_ID, HB_percentile)

# change column names
colnames(testedDS) <- c("lab_ID", "percent_hygienic")

# merge hygieneic behavior back in
ds=merge(x=ds, y=testedDS, all.x = TRUE)



##### Plot HB by varroa load

# Add regression lines
ggplot(ds, aes(x=percent_hygienic, y=varroa_load_mites.100.bees, 
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
ggplot(ds, aes(x=percent_hygienic, y=varroa_load_mites.100.bees)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE) +
  geom_point(size=2) +
  ylab("Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) # size of the text and label ticks




# HB by varroa load model
mod2 <- lm(data=ds, varroa_load_mites.100.bees~percent_hygienic  sampling_event )

summary(mod2)


Anova(mod2)



# plot Hygienic behavior by yard

# take only sampling event 2
dst2 <- ds[ds$sampling_event==2,]




# Hygienic behavior by hive/yard, points number of hives, threshold dotted bar
ggplot(dst2, aes(x=yard, y=percent_hygienic, color=yard)) + 
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
mod3 <- aov(data = dst2, percent_hygienic ~ yard)
summary(mod3)


# print table in increasing order based on some variable
dst2[order(dst2$percent_hygienic, decreasing = TRUE),]





##################################################################################
# selection model
# clean data:

# aggregate mite load by sampling event and yard
modelDF <- ds %>% # operate on the dataframe (ds) and assign to new object 
  group_by(lab_ID, yard) %>% # pick variables to group by
  summarise(
            varroa = mean(varroa_load_mites.100.bees, na.rm=T),
            nosema = mean(nosema_load_spores.bee, na.rm=T),
            hygienic = mean(HB_percentile, na.rm=T),
            honey = mean(honey_removed, na.rm=T),
            weight = mean(october_weight, na.rm=T)
    
            ) 


# impute missing values:
imputedDF <- na_mean(modelDF[,3:7])

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

# create histogram data set of scaled vals
histogram <- data.frame(varroaScaled, nosemaScaled, hygienicScaled, honeyScaled, weightScaled)
histogramLong <- gather(histogram, condition, measurement, varroaScaled, nosemaScaled, hygienicScaled, honeyScaled, weightScaled)


# plot histogram of all scaled values
ggplot(histogramLong ,aes(x=measurement, fill=condition)) + 
  geom_histogram(alpha = 0.5, position = "identity") +
  ylab("Frequency") + # y axis label
  xlab("Parameter Value") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "right") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Parameter") # color pallets option = A-H


# create fitness function
Fitness = (-1 * varroaScaled) + 
  (-1 * nosemaScaled) + 
  (1 * hygienicScaled) + 
  (1 * honeyScaled) + 
  (1 * weightScaled) 

# rescale 0 to 1
imputedDF$Fitness <- range01(Fitness)

# print the sorted data set
imputedDF[order(imputedDF$Fitness, decreasing = TRUE),]


# explicit print command:
print(imputedDF)





plot(sort(imputedDF$Fitness))

hist(imputedDF$Fitness, breaks = 10)


#install.packages("dplyr", "ggplot2", "lme4", "tidyr", "viridis", "car", "imputeTS")

