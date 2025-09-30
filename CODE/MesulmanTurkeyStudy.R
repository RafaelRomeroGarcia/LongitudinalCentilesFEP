
##################################################################################################
# Script for Statistical Analysis and Visualization of Neuroimaging Data
# Author: Claudio Aleman Morillo    email: caleman@us.es
# Date: July 2024
# Version: 1.0

# Introduction:
# This script performs statistical analysis and visualization on neuroimaging data.
# It includes data loading, data transformation, ANOVA, and Tukey's HSD test.
# The script uses generic file paths and column names for ease of adaptation to other datasets.
###################################################################################################

# Description:
#we need to indicate the path to the files with the T values obtained from the mixed linear regressions for 
#the average cognitive functionality and score on the positive (SAPS), negative (SANS) symptom scales
#and on the BPRS scale. These t values are obtained from the execution of the associated scripts. 
#These scripts return their respective csv files, which we have to indicate here. We only have to include the route to them.
#These files are:
#  Global_Cognitive_Functioning_average_T_Values.csv,
#  SAPS_Total_T_Values.csv,
#  SANS_Total_T_Values.csv,
#  BPRS_Total_T_Values.csv.
#From these files we must choose the column to represent, changing the value of the variable 'name_Variable_to_Represent' .
#These columns will be those associated with the centiles depending on the case. 
#For cognition we have: 'Cog_Controls' and 'Cog_Psychosis'.
#For the symptom scales :'SAPS', 'SANS' and 'BPRS' .

#Finally we need to indicate the path of the equivalence file
#between the cortical regions of the Desikan-Killiany atlas and the Mesulam cortical atlas
#in the 'MapDataBase' variable.
#This file is called:
#  'equivalencesMesulan_DK.csv'




###################################################################################################

# Clear the workspace
rm(list=ls())

# Load necessary libraries
library(ggplot2)
library(ggseg)
library(ggseg3d)
library(ggsegExtra)
library(dplyr)
library(lsr)
library(ggpubr)
library(ggsegYeo2011)

# Function to create a bar plot
representation <- function(dataFrame, X, Y, text_title) {
  ggbarplot(dataFrame, x = X, y = Y, fill = "network", add = c("mean_se", "jitter")) + 
    theme(
      axis.line = element_line(),
      axis.text.x = element_text(size = rel(1.8), face = "bold"),
      axis.text.y = element_text(size = rel(1.8), face = "bold"),
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = rel(1.8)),
      axis.title.y = element_text(size = rel(1.8)),
      legend.background = element_rect(fill = "white"),
      legend.key = element_rect(fill = "white"),
      legend.key.size = unit(1.2, "lines"),
      legend.position = "None",
      legend.text = element_text(size = rel(1.8), face = "bold"),
      legend.title = element_text(size = rel(1.8), face = "bold", hjust = 0),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      plot.title = element_text(size = rel(1.2)),
      strip.background = element_rect(fill = "grey90", color = "grey50"),
      strip.text.x = element_text(size = rel(0.8)),
      strip.text.y = element_text(size = rel(0.8), angle = -90)
    ) + 
    ylim(c(-6, 6)) +
    xlab(text_title) +
    ylab('Beta') +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    guides(fill = guide_legend("Areas"))
}

# Function to melt data frame (simplified)
meltDataFrameByColumns <- function(D_dataFrame, id_Var_vectStr) {
  D_dataFrame$namesC <- id_Var_vectStr
  print(D_dataFrame)
  return(D_dataFrame)
}

# Read the mapping data
MapDataBase <- read.csv('path/to/equivalencesMesulan_DK.csv') ##### Modify

# Define column names and file names
columns_Names <- c('dx0.zone', 'dx1.zone', 'zone', 'zone', 'zone')
names(columns_Names) <- c('Cog_Controls', 'Cog_Psychosis', 'SAPS', 'SANS', 'BPRS')
File_Names <- c(
  'path/to/Global_Cognitive_Functioning_average_T_Values.csv',           ##### Modify
  'path/to/Global_Cognitive_Functioning_average_T_Values.csv',           ##### Modify
  'path/to/SAPS_Total_T_Values.csv',                                     ##### Modify
  'path/to/SANS_Total_T_Values.csv',                                     ##### Modify
  'path/to/BPRS_Total_T_Values.csv'                                      ##### Modify
)
names(File_Names) <- c('Cog_Controls', 'Cog_Psychosis', 'SAPS', 'SANS', 'BPRS')

# Variable to represent
name_Variable_to_Represent <- 'SANS'                         ##### Modify values:('Cog_Controls', 'Cog_Psychosis', 'SAPS', 'SANS', 'BPRS')

# Load the appropriate dataset
DataBase <- data.frame(read.csv(File_Names[name_Variable_to_Represent]))
rownames(DataBase) <- DataBase[, 2]
rownames(MapDataBase) <- MapDataBase[,'label1']
DataBase$network <- MapDataBase[rownames(DataBase), 'label2']

# Create the plot
representation(DataBase, 'network', columns_Names[name_Variable_to_Represent], paste0('Loadings ', name_Variable_to_Represent))

# Prepare data for ANOVA
matrixAnovaRaw <- DataBase[, c(columns_Names[name_Variable_to_Represent], 'network')]
colnames(matrixAnovaRaw) <- c('zone', 'network')

# Perform ANOVA
fit <- aov(zone ~ network, data = matrixAnovaRaw)
summary(fit)
etaSquared(fit)
TukeyData <- TukeyHSD(fit)
result <- data.frame(TukeyData$network)

# Print the results
print(result)


