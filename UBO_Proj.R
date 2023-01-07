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
library(cowplot)


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

# create anonymous beekeeper names
ds$anonBeek <- ifelse(ds$beekeeper == "Andrew Munkres", "beekeeper 1",
       ifelse(ds$beekeeper == "Jack Rath", "beekeeper 2", "beekeeper 3"
))

# log transform data
ds$log_june_varroa_load_mites.100.bees <- log10(ds$june_varroa_load_mites.100.bees + 0.0001)
ds$log_august_varroa_load_mites.100.bees <- log10(ds$august_varroa_load_mites.100.bees + 0.0001)
ds$log_june_nosema_load_spores.bee <- log10(ds$june_nosema_load_spores.bee + 1)

# split the data by beekeeper
dsSplit <- split(ds, ds$beekeeper)

## UBO by June Varroa loads continuous 
# Add regression lines
juneBeek <- ggplot(ds, aes(x=assay_score, y=june_varroa_load_mites.100.bees, 
               color=as.character(anonBeek))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2) +
  geom_point(size=3) +
  coord_cartesian(ylim = c(0, 4)) + 
  ylab(NULL) + # y axis label
  xlab(NULL) + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = c(.8,.8)) + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3", "grey"), name="Beekeeper:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))
juneBeek

# discrete analysis
dsNo0 <- ds[!ds$june_varroa_load_mites.100.bees==0,]
boxplot(dsNo0$UBO_binary, dsNo0$june_varroa_load_mites.100.bees)
x = aov(dsNo0$june_varroa_load_mites.100.bees~dsNo0$UBO_binary)
summary(x)


# june varroa for just andrew by ubo
cor.test(dsSplit$`Andrew Munkres`$assay_score, dsSplit$`Andrew Munkres`$june_varroa_load_mites.100.bees, method="spearman", exact = F)
cor.test(dsSplit$`Mike Palmer`$assay_score, dsSplit$`Mike Palmer`$june_varroa_load_mites.100.bees, method="spearman", exact = F)
cor.test(dsSplit$`Jack Rath`$assay_score, dsSplit$`Jack Rath`$june_varroa_load_mites.100.bees, method="spearman", exact = F)


## UBO by August Varroa loads continuous
# Add regression lines
augBeek <- ggplot(ds, aes(x=assay_score, y=august_varroa_load_mites.100.bees, 
               color=as.character(anonBeek))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2) +
  geom_point(size=3) +
  coord_cartesian(ylim = c(0, 4)) + 
  ylab(NULL) + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = c(5, 5)) + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3", "grey"), name="Beekeeper:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))
augBeek

# august varroa load by beekeeper
cor.test(dsSplit$`Andrew Munkres`$assay_score, dsSplit$`Andrew Munkres`$august_varroa_load_mites.100.bees, method="spearman", exact = F)
cor.test(dsSplit$`Mike Palmer`$assay_score, dsSplit$`Mike Palmer`$august_varroa_load_mites.100.bees, method="spearman", exact = F)
cor.test(dsSplit$`Jack Rath`$assay_score, dsSplit$`Jack Rath`$august_varroa_load_mites.100.bees, method="spearman", exact = F)


## UBO June composite continuous 
# Add regression lines
juneVar <- ggplot(ds, aes(x=assay_score, y=june_varroa_load_mites.100.bees)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2, color = "darkturquoise") +
  geom_point(size=3) +
  coord_cartesian(ylim = c(0, 4)) + 
  ylab("June Varroa Load (mites/100 bees)") + # y axis label
  xlab(NULL) + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text("aes(label=lab_ID)", ) +
  guides(color = guide_legend(override.aes = list(label = '')))
juneVar

# june varroa  by ubo
cor.test(ds$assay_score, ds$june_varroa_load_mites.100.bees, method="spearman", exact = F)


## UBO August composite continuous 
# Add regression lines
augVar <- ggplot(ds, aes(x=assay_score, y=august_varroa_load_mites.100.bees)) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2, color = "darkturquoise") +
  geom_point(size=3) +
  coord_cartesian(ylim = c(0, 4)) + 
  ylab("August Varroa Load (mites/100 bees)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_viridis(discrete = TRUE, option="H", name="Time Point:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))
augVar

# june varroa  by ubo
cor.test(ds$assay_score, ds$august_varroa_load_mites.100.bees, method="spearman", exact = F)



# make a multi panel plot
plot_grid(juneVar, juneBeek, augVar, augBeek,
          labels = "AUTO", 
          label_size = 17)





# Hygienic behavior by hive/yard, points number of hives, threshold dotted bar
ggplot(ds, aes(x=anonBeek, y=assay_score, color=beekeeper)) + 
  geom_boxplot(size=1) +
  geom_point(size=3) +
  guides(color = guide_legend(override.aes = list(label = ''))) +
  ylab("Percent Hygienic Behavior") + # y axis label
  xlab("Beekeeper") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "none") + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3", "grey")) +
  geom_hline(yintercept=.6, linetype="dashed",
             color = "black", size=1)



ds %>% # operate on the dataframe (ds_2021) and assign to new object (pltN)
  group_by(beekeeper) %>% # pick variables to group by
  summarise(
    
    mean = mean(UBO_binary, na.rm=T), # mean
    sum = sum(UBO_binary, na.rm=T),
    N = length(UBO_binary), # sample size
    
  ) 



## UBO by June varroa load binary 
boxplot(ds$june_varroa_load_mites.100.bees~ds$UBO_binary)
summary(aov(ds$june_varroa_load_mites.100.bees~ds$UBO_binary))





## UBO by August Varroa loads binary
boxplot(ds$august_varroa_load_mites.100.bees~ds$UBO_binary)
summary(aov(ds$august_varroa_load_mites.100.bees~ds$UBO_binary))






#################################################################################
# Nosema Analysis
#################################################################################


## UBO by June Nosema load continuous 
# Add regression lines
ggplot(ds, aes(x=assay_score, y=june_nosema_load_spores.bee, 
               color=as.character(anonBeek))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2) +
  geom_point(size=3) +
  ylab("June Nosema Load (spores/bee)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3", "grey"), name="Beekeeper:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))

cor.test(dsSplit$`Jack Rath`$assay_score, dsSplit$`Jack Rath`$june_nosema_load_spores.bee, method="spearman", exact = F)
cor.test(dsSplit$`Mike Palmer`$assay_score, dsSplit$`Mike Palmer`$june_nosema_load_spores.bee, method="spearman", exact = F)
cor.test(dsSplit$`Andrew Munkres`$assay_score, dsSplit$`Andrew Munkres`$june_nosema_load_spores.bee, method="spearman", exact = F)



## UBO by June Nosema load continuous wit a log transform
# Add regression lines
ggplot(ds, aes(x=assay_score, y=log_june_nosema_load_spores.bee, 
               color=as.character(anonBeek))) +
  #geom_point(size=0) + 
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 2) +
  geom_point(size=3) +
  ylab("June Nosema Load (spores/bee)") + # y axis label
  xlab("Percent Hygienic Behavior") + # x axis label
  theme_minimal(base_size = 17) + # size of the text and label ticks
  theme(legend.position = "top") + # place the legend at the top
  scale_color_manual(values = c("darkturquoise", "tomato3", "grey"), name="Beekeeper:") + # color pallets option = A-H
  #geom_text(aes(label=lab_ID)) +
  guides(color = guide_legend(override.aes = list(label = '')))


summary(lm(dsSplit$`Andrew Munkres`$log_june_nosema_load_spores.bee ~ dsSplit$`Andrew Munkres`$assay_score))


# remove 0s
dsNos_no0 <- ds[!ds$june_nosema_load_spores.bee==0,] 

# remove NA
dsNos_no0 <- dsNos_no0[!is.na(dsNos_no0$june_nosema_load_spores.bee), ]
dsNos_no0 <- dsNos_no0[!is.na(dsNos_no0$assay_score), ]


#glm with gamma distribution



dsNos_no0$rescaledNosema <- dsNos_no0$june_nosema_load_spores.bee/sum(dsNos_no0$june_nosema_load_spores.bee)


mod <- glmer(data = dsNos_no0, rescaledNosema ~ assay_score + beekeeper + (1 | yard), family = "Gamma", 
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun= 100000)),
              verbose = 1)


Anova(mod)









# june nosema by ubo
cor.test(ds$assay_score, ds$log_june_nosema_load_spores.bee, method="spearman", exact = F)


summary(lm(dsNos_no0$log_june_nosema_load_spores.bee ~ dsNos_no0$assay_score))
plot(dsNos_no0 $assay_score, dsNos_no0 $log_june_nosema_load_spores.bee)

## UBO by June Nosema load binary 
boxplot(dsNos_no0$log_june_nosema_load_spores.bee~dsNos_no0$UBO_binary)
summary(aov(ds$log_june_nosema_load_spores.bee~ds$UBO_binary))

t.test(ds$june_nosema_load_spores.bee~ds$UBO_binary, alternative="greater")

kruskal.test(ds$june_nosema_load_spores.bee~ds$UBO_binary)


 ## UBO by August Nosema load continuous 
# Add regression lines
ggplot(ds, aes(x=assay_score, y=august_nosema_load_spores.bee, 
               color=as.character(beekeeper))) +
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



# here we can create a character variable below 10 is small greater than ten is large


# UBO by frames of brood 
ggplot(ds, aes(x=frames_of_brood, y=assay_score, 
               color=as.character(beekeeper))) +
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






