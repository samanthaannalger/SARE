###########################################################################################
# Functions I Wrote and Use on a Daily Basis
# P. Alexander Burnham
# July 6, 2017

###########################################################################################
# BEGIN FUNCTIONS #########################################################################
###########################################################################################

###########################################################################
# function name: PrelimClean
# description: removes unneed PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=MigVirus){
  
  # take only columns that we want:
  library(dplyr)
  data <- dplyr::select(data, ID, Cq_mean, target_name, quantity_mean, Plate)
  
  # remove duplicate rows
  data<-data[!duplicated(data), ]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="No Sample"),]
  
  # remove Gblock rows from dataframe:
  data<-data[!(data$ID=="Gblock"),]
  
  # remove Gblock rows from dataframe:
  data<-data[!(data$ID=="GBlock"),]
  
  # remove Gblock rows from dataframe:
  data<-data[!(data$ID=="Sample Name"),]
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################



###########################################################################
# function name: VirusNorm
# description: normalizes virus data with dilutions and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 1, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr <- 100
  eluteRNA <- 50
  GITCperbee <- 200
  cDNA_eff <- 0.1
  rxn_vol <- 3
  
  #create column for total_extr_vol
  total_extr_vol <- (GITCperbee * number_bees)
  
  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * data$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  # norm_genomeCopy is 0 if NA
  data$genomeCopy[is.na(data$genomeCopy)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################





###########################################################################
# function name: actinNormal
# description: normalizes virus data with actin values 
# parameters: a data frame with actin values
# returns a dataframe with normalized virus values 
###########################################################################

actinNormal <- function(data=MigVirus){
  
  # pull only actin values out of dataframe
  ActinOnly <- data[which(data$target_name=="ACTIN"),]
  
  # create DF of ACTIN genome copies and lab ID:
  ActinDF <- data.frame(ActinOnly$ID, ActinOnly$Plate, ActinOnly$genomeCopy)
  colnames(ActinDF) <- c("ID", "Plate", "ACT_genomeCopy")
  
  # merge ACTIN dataframe with main dataframe:
  #Need rownames and all.x=TRUE because data frames are different sizes.
  data <- merge(data, ActinDF, by=c("ID", "Plate"), all.x=TRUE)
  
  # find mean of all ACTIN values:
  ActinMean <- mean(ActinOnly$genomeCopy, na.rm = TRUE)
  
  # create column for normalized genome copies per bee:
  data$NormGenomeCopy <- (data$genomeCopy/data$ACT_genomeCopy)*ActinMean
  
  return(data)
}


###########################################################################
# END OF FUNCITON
###########################################################################





###########################################################################
# function name: CT_Threash
# description: creates binary data and makes genome copy 0 if below Ct threash
# parameters: dataframe
###########################################################################

CT_Threash <- function(data=data){
  
  splitDF <- split(data, data$target_name)

  # make DWV norm_genome_copbee 0 if Ct value is > 32.918
  splitDF$DWV$NormGenomeCopy[which(splitDF$DWV$Cq_mean > 32.918)] <- 0
  splitDF$BQCV$NormGenomeCopy[which(splitDF$BQCV$Cq_mean > 32.525)] <- 0
  
  splitDF$DWV$virusBINY  <- ifelse(splitDF$DWV$Cq_mean > 32.918, 0, 1)
  splitDF$BQCV$virusBINY  <- ifelse(splitDF$BQCV$Cq_mean > 32.525, 0, 1)
  
  # merge split dataframe back into "BombSurv" dataframe:
  data <- rbind(splitDF$DWV, splitDF$BQCV, splitDF$IAPV)
  
  # norm_genomeCopy is 0 if NA
  data$virusBINY[is.na(data$virusBINY)] <- 0

return(data)

}

###########################################################################
# END OF FUNCITON
###########################################################################





###########################################################################
# function name: VirusMerger2000
# description: merges binary and viral load data for 3 viruses to data frame
# parameters: needs two data frames (data1 is the main one)
###########################################################################

VirusMerger2000 <- function(data1=DF2, data2=DF2){
  
  MigVirusSplit <- split(data1, data1$target_name)
  
  BQCV <- select(MigVirusSplit$BQCV, sample_name, NormGenomeCopy, virusBINY)
  DWV <- select(MigVirusSplit$DWV, sample_name, NormGenomeCopy, virusBINY)
  IAPV <- select(MigVirusSplit$IAPV, sample_name, NormGenomeCopy, virusBINY)
  
  sample_name <- BQCV$sample_name
  BQCVload <- BQCV$NormGenomeCopy
  BQCVbinary <- BQCV$virusBINY
  DWVload <- DWV$NormGenomeCopy
  DWVbinary <- DWV$virusBINY
  IAPVload <- IAPV$NormGenomeCopy
  IAPVbinary <-  IAPV$virusBINY
  
  DF <- data.frame(sample_name, BQCVload, DWVload, IAPVload, DWVbinary, IAPVbinary, BQCVbinary)
  
  data <- merge(data2, DF, by = c("sample_name"), all.x = T)
  
  return(data)
}

###########################################################################
# END OF FUNCITON
###########################################################################

